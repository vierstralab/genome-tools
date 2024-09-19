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


class base_extractor(object):
    def __init__(self, filename, **kwargs):
        self.filename = filename

    def __getitem__(self, i):
        raise NotImplementedError


# ------------------------


class fasta_extractor(base_extractor):
    def __init__(self, filename, **kwargs):
        super(fasta_extractor, self).__init__(filename, **kwargs)

        self.fasta = pysam.FastaFile(filename)

    def __getitem__(self, interval):
        seq = self.fasta.fetch(interval.chrom, interval.start, interval.end)
        # implement any post processing steps here
        return seq

    def __del__(self):
        self.close()

    def close(self):
        if self.fasta and self.fasta.is_open():
            self.fasta.close()


# ------------------------


class tabix_iter(object):
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


class tabix_extractor(base_extractor):
    def __init__(self, filename, header_char="#", columns=None, **kwargs):
        """ """
        super(tabix_extractor, self).__init__(filename, **kwargs)

        self.tabix = pysam.TabixFile(filename)


        with gzip.open(filename, "rt") as f:
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
        ret = pd.read_table(
            tabix_iter(self.tabix, interval), header=None, index_col=None
        )
        ret.columns = self.columns
        return ret

    def __del__(self):
        self.close()

    def close(self):
        if getattr(self, "tabix", None) and self.tabix.is_open():
            self.tabix.close()


# ------------------------


class bigwig_extractor(base_extractor):
    def __init__(self, filename, **kwargs):
        """ """
        super(bigwig_extractor, self).__init__(filename, **kwargs)

        self.bw = pbw.open(filename)

    def __getitem__(self, interval):
        out = self.bw.values(interval.chrom, interval.start, interval.end, numpy=True)
        return out

    def __del__(self):
        self.close()

    def close(self):
        if getattr(self, "bw", None):
            self.bw.close()


# ------------------------


class d4_extractor(base_extractor):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

    def __getitem__(self, i):
        raise NotImplementedError

    def __del__(self):
        self.close()

    def close(self):
        pass
