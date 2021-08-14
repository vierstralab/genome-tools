"""Classes to extract genomics data from different file formats

Note that many of these filehandlers are not thread safe,
so they must be opened separately on different threads when
using multiprocessing.
"""
import numpy as np
import pysam
import pyBigWig as pbg

class base_extractor(object):
    def __init__(self, filename, **kwargs):
        self.filename = filename

    def __getitem__(self, i):
        raise NotImplementedError 

class fasta_extractor(base_extractor):
    def __init__(self, filename, *kwargs):
        """
        """
        super(fasta_extractor, self).__init__(filename, **kwargs)
        
        self.fasta = pysam.FastaFile(filename)

    def __getitem__(self, interval):
        """
        """
        seq = self.fasta.fetch(interval.chrom, interval.start, interval.end)
        # implement any post processing steps here
        return seq

    def __del__(self):
        self.close()

    def close(self):
        if self.fasta and self.fasta.is_open():
            self.fasta.close()

class tabix_extractor(base_extractor):
    def __init__(self, filename, col=4, dtype=float, **kwargs):
        """
        """
        super(tabix_extractor, self).__init__(filename, **kwargs)

        self.tabix = pysam.TabixFile(filename)
        self.col = col
        self.dtype = dtype

    def __getitem__(self, interval):
        out = np.zeros(len(interval), dtype=self.dtype)
        for e in self.tabix.fetch(interval.chrom, interval.start, 
            interval.end, parser=pysam.asTuple()):
            idx = int(e[1])-interval.start
            out[idx] = self.dtype(e[self.col])
        return out

    def __del__(self):
        self.close()

    def close(self):
        if self.tabix and self.tabix.is_open():
            self.tabix.close()
    
class bigwig_extractor(base_extractor):
    def __init__(self, filename, **kwargs):
        """
        """
        super(bigwig_extractor, self).__init__(filename, **kwargs)

        self.bw = pbw.open(filename)

    def __getitem__(self, interval):
        out = self.bw.values(interval.chrom, interval.start, interval.end)
        return out

    def __del__(self):
        self.close()

    def close(self):
        if self.bw:
            self.bw.close()