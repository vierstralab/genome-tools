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
import numpy as np

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
    """Wraps tabix fetch to return an iterator that can be used with pandas"""

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
                'chrom': record.chrom,
                'pos': record.pos,
                'id': record.id,
                'ref': record.ref,
                'alt': alt,
                'gt': [y['GT'] for y in record.samples.values()],
                'indiv_id': list(record.samples.keys()),
            }
            variants.append(variant)
        df = pd.DataFrame(variants).explode(['gt', 'indiv_id'])
        return df

    def close(self):
        if getattr(self, "variant_vcf", None) and self.variant_vcf.is_open():
            self.variant_vcf.close()


class OverlappingReadsExtractor(BaseExtractor):
    def __init__(self, filename, **kwargs):
        """
        Extracts data from a cram/bam file.

        Parameters
        ----------
        filename : str
            Path to the cram/bam file.
        """
        super().__init__(filename, **kwargs)

        self.cram = pysam.AlignmentFile(filename)

    def __getitem__(self, interval):
        reads_1 = {}
        reads_2 = {}

        # Go into BAM file and get the reads
        for pileupcolumn in self.cram.pileup(interval.chrom, interval.start, interval.end, maxdepth=10000, truncate=True, stepper="nofilter"):
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del or pileupread.is_refskip:
                    print('refskip or del ', pileupread.alignment.query_name, file=sys.stderr)
                    continue

                if pileupread.alignment.is_read1:
                    reads_1[pileupread.alignment.query_name] = pileupread
                else:
                    reads_2[pileupread.alignment.query_name] = pileupread

        # All reads that overlap the region; unique set
        read_pairs = set(reads_1.keys()) | set(reads_2.keys())

        return reads_1, reads_2, read_pairs


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
            line = f.readline().strip('\n')
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
                TabixIter(self.tabix, interval), header=None, index_col=None, names=self.columns, **self.kwargs
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
            filters=[('chrom', '==', interval.chrom)],
            engine='pyarrow',
            columns=self.columns,
        ).iloc[interval.start:interval.end]

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
