# Copyright 2015 Jeff Vierstra

import numpy as np


class genomic_interval(object):
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

    def __str__(self):
        """Returns a string-formated version of the element
        for printing
        """
        return "\t".join([str(x) for x in [self.chrom, self.start, self.end]])

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
        return genomic_interval(
            self.chrom, self.start - x, self.end + x, self.name, self.strand
        )

    def shift(self, x, inplace=False):
        """Shifts the coordinates"""
        if inplace:
            self.start += x
            self.end += x
            return self
        return genomic_interval(
            self.chrom, self.start + x, self.end + x, self.name, self.strand
        )


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
    result = [genomic_interval(df_row['chrom'], df_row['start'], df_row['end'], **{col: df_row[col] for col in extra_columns}) for _, df_row in df.iterrows()]
    return result
