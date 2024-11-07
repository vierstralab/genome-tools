# Copyright 2015 Jeff Vierstra

import numpy as np


class GenomicInterval(object):
    """Class that implements BED-style object"""

    def __init__(
        self, chrom, start, end, name=".", **kwargs
    ):
        self.chrom = str(chrom)
        self.start = int(start)
        self.end = int(end)
        self.name = str(name)
        self._extra_fields = []
        for key, value in kwargs.items():
            self._extra_fields.append(key)
            setattr(self, key, value)

    def from_ucsc(ucsc_str):
        """Parses a UCSC-style string to a GenomicInterval object"""
        chrom, rest = ucsc_str.split(":")
        start, end = rest.split("-")
        return GenomicInterval(chrom, int(start.replace(',', '')), int(end.replace(',', '')))

    def __len__(self):
        """Length of element"""
        return self.end - self.start

    def __repr__(self):
        """Returns a string-formated version of the element
        for printing
        """
        extra_fields = f", extra_fields={tuple(self._extra_fields)}" if len(self._extra_fields) > 0 else ""
        return f'{self.__class__.__name__}({self.to_ucsc()}{extra_fields})'

    def to_ucsc(self):
        """Prints out a UCSC version of interval. Note
        that the start coordinate is shifted by 1 to
        reflect 1-based"""
        return f"{self.chrom}:{self.start}-{self.end}"
    
    @property
    def _extra_kwargs(self):
        return {key: getattr(self, key) for key in self._extra_fields}

    def widen(self, x, inplace=False):
        """Expands the coordinates"""
        new_start = max(self.start - x, 0)
        new_end = self.end + x
        if new_end < new_start:
            mid = (new_start + new_end) // 2
            new_start = mid
            new_end = mid + 1
        if inplace:
            self.start = new_start
            self.end = new_end
            return self
        return GenomicInterval(
            self.chrom, new_start, new_end, self.name, **self._extra_kwargs
        )
    
    def zoom(self, zoom_factor, inplace=False):
        """
        Zoom in/out on the coordinates
        Positive zoom is interpreted as zooming in
        Negative zoom is interpreted as zooming out
        """
        assert zoom_factor != 0, "Zoom factor cannot be zero"
        if zoom_factor < 0:
            zoom_factor = -1 / zoom_factor
        
        padding = int(len(self) * ((1 / zoom_factor) - 1) / 2)
        return self.widen(padding, inplace=inplace)

    def shift(self, x, inplace=False):
        """Shifts the coordinates"""
        if inplace:
            self.start += x
            self.end += x
            return self
        return GenomicInterval(
            self.chrom, self.start + x, self.end + x, self.name, **self._extra_kwargs
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
    
    def zoom(self, zoom_factor, inplace=False):
        raise NotImplementedError("Cannot zoom a variant interval")


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
