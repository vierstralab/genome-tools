"""Classes to extract genomics data from different file formats

Note that many of these filehandlers are not thread safe,
so they must be opened separately on different threads when
using multiprocessing.
"""
import numpy as np
import pysam
import pyBigWig as pbw
import pandas as pd
import gzip



class BaseExtractor:
    def __init__(self, filename, **kwargs):
        self.filename = filename

    def __getitem__(self, i):
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

    def __init__(self, tabix, interval):
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
            Additional arguments to pass to the base class.
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
                TabixIter(self.tabix, interval), header=None, index_col=None
            )
        except pd.errors.EmptyDataError:
            ret = pd.DataFrame(columns=self.columns)
        ret.columns = self.columns
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
