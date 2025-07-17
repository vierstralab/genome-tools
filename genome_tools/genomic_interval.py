# Copyright 2015 Jeff Vierstra

from typing import List
import pandas as pd
from tqdm import tqdm


class GenomicInterval(object):
    """Class that implements BED-style object"""

    def __init__(self, chrom, start, end, name=".", **kwargs):
        super().__setattr__('_initialized', False)
        self.chrom = str(chrom)
        self.start = int(start)
        self.end = int(end)
        self.name = str(name)
        self.req_fields = ["chrom", "start", "end", "name"]
        self.extra_fields = []
        for key, value in kwargs.items():
            self.extra_fields.append(key)
            setattr(self, key, value)
        super().__setattr__('_initialized', True)

    def from_ucsc(ucsc_str):
        """Parses a UCSC-style string to a GenomicInterval object"""
        chrom, rest = ucsc_str.split(":")
        start, end = rest.split("-")
        return GenomicInterval(chrom, int(start.replace(',', '')), int(end.replace(',', '')))

    def __len__(self):
        """Length of element"""
        return self.end - self.start

    def __setattr__(self, name, value):
        if getattr(self, '_initialized', False):
            reserved_fields = ('extra_fields', 'req_fields')
            if name in reserved_fields:
                raise AttributeError(f"Cannot set reserved fields directly {reserved_fields}")
            if not hasattr(self, name):
                self.extra_fields.append(name)
        super().__setattr__(name, value)

    def __repr__(self):
        """Returns a string-formated version of the element
        for printing
        """
        extra_fields = f", extra_fields={tuple(self.extra_fields)}" if len(self.extra_fields) > 0 else ""
        return f'{self.__class__.__name__}({self.to_ucsc()}{extra_fields})'

    def to_ucsc(self):
        """Prints out a UCSC version of interval. Note
        that the start coordinate is shifted by 1 to
        reflect 1-based"""
        return f"{self.chrom}:{self.start}-{self.end}"
    
    @property
    def extra_kwargs(self):
        return {key: getattr(self, key) for key in self.extra_fields}
    
    @property
    def center(self):
        return self.widen(-len(self))
    
    def __add__(self, x: int):
        return self.shift(x)
    
    def __iadd__(self, x):
        return self.shift(x, inplace=True)
    
    def __sub__(self, x: int):
        return self.shift(-x)
    
    def __isub__(self, x: int):
        return self.shift(-x, inplace=True)
    
    def __mul__(self, x: float):
        return self.zoom(x)
    
    def __imul__(self, x: float):
        return self.zoom(x, inplace=True)
    
    def __truediv__(self, x: float):
        return self.zoom(1 / x)
    
    def __itruediv__(self, x: float):
        return self.zoom(1 / x, inplace=True)

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
            self.chrom, new_start, new_end, self.name, **self.extra_kwargs
        )
    
    def zoom(self, zoom_factor, inplace=False):
        """
        Zoom in/out on the coordinates
        Positive zoom is interpreted as zooming in
        Negative zoom is interpreted as zooming out
        """
        if not isinstance(zoom_factor, (int, float)):
            raise TypeError("zoom factor must be an int or float")
        assert zoom_factor != 0, "Zoom factor cannot be zero"
        if zoom_factor < 0:
            zoom_factor = -1 / zoom_factor
        
        padding = round(len(self) * ((1 / zoom_factor) - 1) / 2)
        return self.widen(padding, inplace=inplace)

    def shift(self, x: int, inplace=False):
        if not isinstance(x, int):
            raise TypeError("shift must be an integer")
        """Shifts the coordinates"""
        if inplace:
            self.start += x
            self.end += x
            return self
        return GenomicInterval(
            self.chrom, self.start + x, self.end + x, self.name, **self.extra_kwargs
        )
    

class VariantInterval(GenomicInterval):
    def __init__(self, chrom, start, end, ref=None, alt=None, value=None, **kwargs):
        assert end - start == 1
        super().__init__(chrom, start, end, **kwargs)
        super().__setattr__('_initialized', False)
        self.ref = ref
        self.alt = alt
        self.value = value
        self.req_fields += ["ref", "alt", "value"]
        self.pos = end
        super().__setattr__('_initialized', True)

    def to_ucsc(self):
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"
    
    def widen(self, x, inplace=False):
        raise NotImplementedError("Cannot widen a variant interval")
    
    def shift(self, x, inplace=False):
        raise NotImplementedError("Cannot shift a variant interval")
    
    def zoom(self, zoom_factor, inplace=False):
        raise NotImplementedError("Cannot zoom a variant interval")


def filter_df_to_interval(df: pd.DataFrame, interval: GenomicInterval):
    chromosome_col = 'chrom' if 'chrom' in df.columns else '#chr'
    return df.loc[df[chromosome_col] == interval.chrom].query(
        f'end >= {interval.start} & start < {interval.end}'
    )


def _sanitize_df(df: pd.DataFrame, interval: GenomicInterval):
    df = df.rename(
        columns={'#chr': 'chrom'}
    )
    if interval is not None:
        df = filter_df_to_interval(df, interval)
    return df


def df_to_genomic_intervals(df: pd.DataFrame, interval: GenomicInterval = None, extra_columns=(), verbose=False):
    df = _sanitize_df(df, interval)
    iterator = df.itertuples(index=False, name=None)
    if verbose:
        iterator = tqdm(iterator, total=len(df), desc="Converting df to intervals")
    get_col_idx = {col: i for i, col in enumerate(df.columns)}
    result = [
        GenomicInterval(
            row[get_col_idx['chrom']],
            row[get_col_idx['start']],
            row[get_col_idx['end']],
            **{col: row[get_col_idx[col]] for col in extra_columns}
        )
        for row in iterator
    ]
    return result


def df_to_variant_intervals(df: pd.DataFrame, interval: GenomicInterval=None, extra_columns=()):
    df = _sanitize_df(df, interval)
    result = [
        VariantInterval(
            df_row['chrom'],
            df_row['start'],
            df_row['end'],
            ref=df_row['ref'],
            alt=df_row['alt'],
            **{col: df_row[col] for col in extra_columns}
        ) for _, df_row in df.iterrows()
    ]
    return result


def genomic_intervals_to_df(genomic_intervals: List[GenomicInterval]):
    assert len(genomic_intervals) > 0, "No intervals to convert"
    req_fields = genomic_intervals[0].req_fields
    extra_fields = genomic_intervals[0].extra_fields
    assert all(interval.extra_fields == extra_fields for interval in genomic_intervals), "Not all _extra_fields of genomic_intervals are identical!"
    fields = req_fields + extra_fields
    return pd.DataFrame(
        [
            [getattr(interval, field) for field in fields]
            for interval in genomic_intervals
        ],
        columns=fields
    )
