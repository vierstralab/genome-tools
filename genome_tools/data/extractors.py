"""Classes to extract genomics data from different file formats

Note that many of these filehandlers are not thread safe,
so they must be opened separately on different threads when
using multiprocessing.
"""

import pysam
import pyBigWig as pbw
import pandas as pd
import gzip
import sys

from genome_tools import GenomicInterval, VariantInterval


class BaseExtractor:
    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.kwargs = kwargs

    def __getitem__(self, interval: GenomicInterval):
        raise NotImplementedError

    def close(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __del__(self):
        self.close()


# ------------------------


class FastaExtractor(BaseExtractor):
    def __init__(self, filename, **kwargs):
        super(FastaExtractor, self).__init__(filename, **kwargs)

        self.fasta = pysam.FastaFile(filename)

    def __getitem__(self, interval):
        seq = self.fasta.fetch(interval.chrom, interval.start, interval.end)
        # implement any post processing steps here
        return seq

    def close(self):
        if self.fasta is not None and not self.fasta.closed:
            self.fasta.close()


# ------------------------


class TabixIter(object):
    """Wraps tabix fetch to return an iterator that can be used with Pandas"""

    def __init__(self, tabix: pysam.TabixFile, interval: GenomicInterval):
        self.iter = tabix.fetch(interval.chrom, interval.start, interval.end)

    def read(self, n=0):
        try:
            return next(self.iter) + "\n"
        except StopIteration:
            return ""

    def __iter__(self):
        return self

    def __next__(self):
        return self.read()


class VariantGenotypeExtractor(BaseExtractor):
    def __init__(self, filename, **kwargs):
        """
        Extracts data from a VCF file. The file must be bgzip compressed and indexed with tabix.

        Parameters
        ----------
        filename : str
            Path to the VCF file.
        **kwargs : dict
            Additional arguments passed to `pd.read_table`.
        """
        super().__init__(filename, **kwargs)

        self.variant_vcf = pysam.VariantFile(filename)

    def __getitem__(self, interval):
        records = self.variant_vcf.fetch(interval.chrom, interval.start, interval.end)
        variants = []

        for record in records:
            assert len(record.alts) == 1
            alt = record.alts[0]
            if isinstance(interval, VariantInterval):
                alt_allele = interval.alt
                assert record.ref == interval.ref
                if alt != alt_allele:
                    continue

            variant = {
                "chrom": record.chrom,
                "pos": record.pos,
                "id": record.id,
                "ref": record.ref,
                "alt": alt,
                "gt": [y["GT"] for y in record.samples.values()],
                "indiv_id": list(record.samples.keys()),
            }
            variants.append(variant)
        df = pd.DataFrame(variants).explode(["gt", "indiv_id"])
        return df

    def close(self):
        if getattr(self, "variant_vcf", None) and self.variant_vcf.is_open():
            self.variant_vcf.close()


# ------------------------


# TODO: eliminate double for loop over reads in _get_all_pileup_reads and _aggregate_into_pairs
class AllelicReadsExtractor(BaseExtractor):
    def __init__(self, filename, variant_read_end_offset=3, **kwargs):
        """
        Extracts data from a cram/bam file. Keep in mind that no adjustments are made for cutcounts offsets (0, -1) for DNase I. Use FootprintTools BamFileExtractor if you need that.

        Parameters
        ----------
        filename : str
            Path to the cram/bam file.
        variant_read_end_offset : int
            Minimum distance from read end to variant position to consider read unbiased.
        **kwargs : dict
        """
        super().__init__(filename, **kwargs)

        self.variant_read_end_offset = variant_read_end_offset
        self.cram = pysam.AlignmentFile(filename, **kwargs)

    @staticmethod
    def get_5p_offset(pileupread: pysam.PileupRead):
        """
        Returns position of variant relative to 5' of read
        """
        if pileupread.query_position is None:  # pileup overlaps deletion
            return None
        elif pileupread.alignment.is_reverse:
            return pileupread.alignment.query_length - pileupread.query_position
        else:
            return pileupread.query_position + 1

    def check_bias(self, pileupread):
        # if get_base_quality(pileupread)<baseq:
        # 	raise ReadBiasError()
        return self.get_5p_offset(pileupread) > self.variant_read_end_offset

    def _validate_read(
        self, pileupread: pysam.PileupRead, variant_interval: VariantInterval
    ):
        if pileupread.is_del or pileupread.is_refskip:
            return pileupread, "N", "deletion_or_refskip"
        overlaps_variant = pileupread.query_position is not None
        if not overlaps_variant:
            return pileupread, "N", "no_overlap"

        if not self.check_bias(pileupread):
            return pileupread, "N", "bias"

        read_allele = pileupread.alignment.query_sequence[pileupread.query_position]
        if read_allele not in (variant_interval.ref, variant_interval.alt):
            return pileupread, "N", "genotype_error"
        return pileupread, read_allele, "valid"

    def _get_all_pileup_reads(self, variant_interval: VariantInterval):
        assert isinstance(variant_interval, VariantInterval), (
            "variant_interval must be an instance of genome_tools.VariantInterval"
        )

        # Go into BAM file and get the reads
        for pileupcolumn in self.cram.pileup(
            variant_interval.chrom,
            variant_interval.start,
            variant_interval.end,
            maxdepth=10000,
            truncate=True,
            stepper="nofilter",
        ):
            pileupcolumn: pysam.PileupColumn
            for pileupread in pileupcolumn.pileups:
                pileupread: pysam.PileupRead
                yield self._validate_read(pileupread, variant_interval)

    @staticmethod
    def _aggregate_into_pairs(reads: dict):
        unique_reads = {}
        for pileupread, read_allele, status in reads:
            pileupread: pysam.PileupRead
            read_name = pileupread.alignment.query_name
            if read_name not in unique_reads:
                unique_reads[read_name] = {
                    "status": status,
                    "read_allele": read_allele,
                    "pileupread": pileupread,
                }
            else:  # goes to this branch if mate exists
                # check that both reads have compatible status.
                mate_status = unique_reads[read_name]["status"]
                if mate_status == status:
                    if status == "valid":
                        mate_allele = unique_reads[read_name]["read_allele"]
                        if mate_allele != read_allele:
                            unique_reads[read_name]["status"] = "different_genotypes"
                            unique_reads[read_name]["read_allele"] = (
                                f"{mate_allele};{read_allele}"
                            )

                elif status == "valid":
                    unique_reads[read_name]["status"] = status
                elif mate_status == "valid":
                    continue
                else:
                    unique_reads[read_name]["status"] = f"{mate_status};{status}"

        return unique_reads

    def extract_allelic_reads(self, variant_interval: VariantInterval):
        assert isinstance(variant_interval, VariantInterval), (
            "variant_interval must be an instance of genome_tools.VariantInterval"
        )
        reads = self._get_all_pileup_reads(variant_interval)
        unique_reads = self._aggregate_into_pairs(reads)
        return unique_reads

    @staticmethod
    def genomic_interval_from_read(pileupread: pysam.PileupRead, base):
        read = pileupread.alignment
        return GenomicInterval(
            chrom=read.reference_name,
            start=read.reference_start,
            end=read.reference_end,
            is_reverse=read.is_reverse,
            base=base,
        )

    def __getitem__(self, variant_interval: VariantInterval):
        allelic_reads = self.extract_allelic_reads(variant_interval)
        return [
            self.genomic_interval_from_read(v["pileupread"], v["read_allele"])
            for v in allelic_reads.values()
            if v["status"] == "valid"
        ]

    def close(self):
        if getattr(self, "cram", None) and not self.cram.closed:
            self.cram.close()


class TabixExtractor(BaseExtractor):
    def __init__(self, filename, header_char="#", columns=None, skiprows=0, **kwargs):
        """
        Extracts data from a tabix file. The file must be bgzip compressed and indexed with tabix.

        Parameters
        ----------
        filename : str
            Path to the tabix file.
        header_char : str
            Character that indicates the start of a header line. Default is '#'.
        columns : list
            List of column names to use. If None, the first line of the file is used as the header.
        skiprows : int
            Number of rows to skip at the beginning of the file. Default is 0.
        **kwargs : dict
            Additional arguments passed to `pd.read_table`.
        """
        super(TabixExtractor, self).__init__(filename, **kwargs)

        self.tabix = pysam.TabixFile(filename)

        with gzip.open(filename, "rt") as f:
            for _ in range(skiprows):
                next(f)
            line = f.readline().strip("\n")
            if columns is None:
                if line.startswith(header_char):
                    self.columns = line.split("\t")
                else:
                    self.columns = [i for i in range(len(line.split("\t")))]
            else:
                assert len(columns) == len(line.split("\t"))
                self.columns = columns

    def __getitem__(self, interval):
        try:
            ret = pd.read_table(
                TabixIter(self.tabix, interval),
                header=None,
                index_col=None,
                names=self.columns,
                **self.kwargs,
            )
        except pd.errors.EmptyDataError:
            ret = pd.DataFrame(columns=self.columns)

        return ret

    def close(self):
        if getattr(self, "tabix", None) and self.tabix.is_open():
            self.tabix.close()


# ------------------------


class ChromParquetExtractor(BaseExtractor):
    def __init__(self, filename, columns=None, **kwargs):
        """
        Extracts data from a single bp resolution parquet file partitioned by chromosome
        """
        super(ChromParquetExtractor, self).__init__(filename, **kwargs)
        self.columns = columns

    def __getitem__(self, interval):
        return pd.read_parquet(
            self.filename,
            filters=[("chrom", "==", interval.chrom)],
            engine="pyarrow",
            columns=self.columns,
        ).iloc[interval.start : interval.end]

    def close(self):
        pass


# ------------------------


class BigwigExtractor(BaseExtractor):
    def __init__(self, filename, **kwargs):
        """ """
        super(BigwigExtractor, self).__init__(filename, **kwargs)

        self.bw = pbw.open(filename)

    def __getitem__(self, interval):
        out = self.bw.values(interval.chrom, interval.start, interval.end, numpy=True)
        return out

    def close(self):
        if getattr(self, "bw", None):
            self.bw.close()


# ------------------------


class D4Extractor(BaseExtractor):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

    def __getitem__(self, i):
        raise NotImplementedError

    def close(self):
        pass
