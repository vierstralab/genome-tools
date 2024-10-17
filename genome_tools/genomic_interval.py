# Copyright 2015 Jeff Vierstra

import numpy as np


class GenomicInterval(object):
    """Class that implements BED-style object"""

    def __init__(
        self, chrom, start, end, name=".", score=None, strand=None, extra=None, **kwargs
    ):
        self.chrom = str(chrom)
        self.start = int(start)
        self.end = int(end)
        self.name = str(name)
        self.score = score
        self.strand = strand
        self.extra = extra
        for key, value in kwargs.items():
            setattr(self, key, value)

    def __len__(self):
        """Length of element"""
        return self.end - self.start

    def __repr__(self):
        """Returns a string-formated version of the element
        for printing
        """
        return f'{self.__class__.__name__}({self.to_ucsc()})'

    def to_ucsc(self):
        """Prints out a UCSC version of interval. Note
        that the start coordinate is shifted by 1 to
        reflect 1-based"""
        return f"{self.chrom}:{self.start}-{self.end}"

    def widen(self, x, inplace=False):
        """Expands the coordinates"""
        if inplace:
            self.start -= x
            self.end += x
            return self
        return GenomicInterval(
            self.chrom, self.start - x, self.end + x, self.name, self.strand
        )

    def shift(self, x, inplace=False):
        """Shifts the coordinates"""
        if inplace:
            self.start += x
            self.end += x
            return self
        return GenomicInterval(
            self.chrom, self.start + x, self.end + x, self.name, self.strand
        )
    

class VariantInterval(GenomicInterval):
    def __init__(self, chrom, start, end, ref=None, alt=None, value=None, **kwargs):
        assert end - start == 1
        super().__init__(chrom, start, end, **kwargs)
        self.ref = ref
        self.alt = alt
        self.pos = end
        self.value = value

    def to_ucsc(self):
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"
    
    def widen(self, x, inplace=False):
        raise NotImplementedError("Cannot widen a variant interval")
    
    def shift(self, x, inplace=False):
        raise NotImplementedError("Cannot shift a variant interval")


def filter_df_to_interval(df, interval):
    chromosome_col = 'chrom' if 'chrom' in df.columns else '#chr'
    return df.loc[df[chromosome_col] == interval.chrom].query(
        f'end >= {interval.start} & start < {interval.end}'
    )



def df_to_genomic_intervals(df, interval=None, extra_columns=()):
    df = df.rename(
        columns={'#chr': 'chrom'}
    )
    if interval is not None:
        df = filter_df_to_interval(df, interval)
    result = [GenomicInterval(df_row['chrom'], df_row['start'], df_row['end'], **{col: df_row[col] for col in extra_columns}) for _, df_row in df.iterrows()]
    return result


def df_to_variant_intervals(df, interval=None, extra_columns=()):
    df = df.rename(
        columns={'#chr': 'chrom'}
    )
    if interval is not None:
        df = filter_df_to_interval(df, interval)
    result = [VariantInterval(df_row['chrom'], df_row['start'], df_row['end'], ref=df_row['ref'], alt=df_row['alt'], **{col: df_row[col] for col in extra_columns}) for _, df_row in df.iterrows()]
    return result