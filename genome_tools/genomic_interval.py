# Copyright 2015 Jeff Vierstra

from typing import List
import pandas as pd
from tqdm import tqdm
import warnings


class GenomicInterval:
    """Class that implements BED-style object"""
    required_fields = ["chrom", "start", "end"]

    def __init__(self, chrom, start, end, name=".", **kwargs):
        super().__setattr__('_initialized', False)
        self.chrom = str(chrom)
        self.start = int(start)
        self.end = int(end)
        self.name = str(name)
        self.extra_fields = []
        for key, value in kwargs.items():
            self.extra_fields.append(key)
            setattr(self, key, value)
        self.extra_fields.append("name")
        super().__setattr__('_initialized', True)

    @staticmethod
    def from_ucsc(ucsc_str):
        """Parses a UCSC-style string to a GenomicInterval object"""
        chrom, rest = ucsc_str.split(":")
        start, end = rest.split("-")
        return GenomicInterval(chrom, int(start.replace(',', '')), int(end.replace(',', '')))

    def copy(self):
        """Returns a copy of the object"""
        return self.__class__(**self._all_kwargs)

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
    
    def _get_kwargs(self, keys):
        return {
            key: getattr(self, key)
            for key in keys
        }
    
    @property
    def _all_kwargs(self):
        return {
            **self._get_kwargs(self.required_fields + self.extra_fields),
        }
    
    @property
    def extra_kwargs(self):
        return self._get_kwargs(self.extra_fields)
    
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

    def widen(self, x=None, *, inplace=False, left=None, right=None):
        """Expands the coordinates"""
        if left is not None:
            assert right is not None, "If left is specified, right must also be specified"
            if x is not None:
                raise ValueError("Either x or both left and right must be specified, not both")
        else:
            left = right = x

        new_start = max(self.start - left, 0)
        new_end = self.end + right
        if new_end < new_start:
            mid = (new_start + new_end) // 2
            new_start = mid
            new_end = mid + 1
        if inplace:
            self.start = new_start
            self.end = new_end
            return self
        return GenomicInterval(
            self.chrom,
            new_start,
            new_end,
            **self.extra_kwargs
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
        return self.widen(left=-x, right=x, inplace=inplace)
    
    def to_variant_interval(self, ref=None, alt=None):
        kwargs = {"ref": ref, "alt": alt}
        assert len(self) == 1, "Can only convert to VariantInterval if length is 1"

        for key, val in kwargs.items():
            if val is None:
                if hasattr(self, key):
                    kwargs[key] = getattr(self, key)
                else:
                    raise ValueError(f"Cannot convert to VariantInterval without {key} specified")

        return VariantInterval(
            **self._all_kwargs,
            ref=kwargs["ref"],
            alt=kwargs["alt"],
        )
    
    def to_genomic_interval(self):
        """Returns self"""
        return self.copy()
    
    def vi(self, ref=None, alt=None):
        """Shortcut for to_variant_interval"""
        return self.to_variant_interval(ref=ref, alt=alt)
    
    def gi(self):
        """Shortcut for to_genomic_interval"""
        return self.to_genomic_interval()

    def overlaps(self, other: 'GenomicInterval'):
        """Returns whether two intervals overlap"""
        assert isinstance(other, GenomicInterval), "other must be a GenomicInterval"
        if self.chrom != other.chrom:
            return False
        return self.end > other.start and self.start < other.end
    

class VariantInterval(GenomicInterval):
    required_fields = GenomicInterval.required_fields + ["ref", "alt"]

    def __init__(self, chrom, start, end, ref, alt, **kwargs):
        assert end - start == 1
        super().__init__(chrom, start, end, **kwargs)
        super().__setattr__('_initialized', False)
        self.ref = str(ref)
        self.alt = str(alt)
        self.pos = end
        super().__setattr__('_initialized', True)

    def to_str(self):
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"

    @staticmethod
    def from_str(variant_str: str):
        chrom, pos, ref, alt = variant_str.split(":")
        return VariantInterval(chrom, int(pos) - 1, int(pos), ref, alt)
    
    def to_genomic_interval(self):
        return GenomicInterval(**self._all_kwargs)
    
    def to_variant_interval(self):
        """Returns self"""
        return self
    
    def widen(self, x, inplace=False):
        raise NotImplementedError("Cannot widen a variant interval. Convert to GenomicInterval with gi() method first.")
    
    def shift(self, x, inplace=False):
        raise NotImplementedError("Cannot shift a variant interval. Convert to GenomicInterval with gi() method first.")
    
    def zoom(self, zoom_factor, inplace=False):
        raise NotImplementedError("Cannot zoom a variant interval. Convert to GenomicInterval with gi() method first.")


class genomic_interval(GenomicInterval):

    _warned = False
    def __init__(self, chrom, start, end, name=".", **kwargs):
        if not genomic_interval._warned:
            warnings.warn(
                "genomic_interval is deprecated; use GenomicInterval instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            genomic_interval._warned = True
        super().__init__(chrom, start, end, name, **kwargs)


def filter_df_to_interval(df: pd.DataFrame, interval: GenomicInterval, strict: bool = False):
    chromosome_col = 'chrom' if 'chrom' in df.columns else '#chr'
    df_slice = df.loc[df[chromosome_col] == interval.chrom].query(
        f'end > {interval.start} & start < {interval.end}'
    )
    if strict:
        df_slice = df_slice.query(
            f'start >= {interval.start} & end <= {interval.end}'
        )
    if isinstance(interval, VariantInterval):
        if 'ref' in df.columns and 'alt' in df.columns:
            df_slice = df_slice.query(f'ref == "{interval.ref}" & alt == "{interval.alt}"')
        else:
            raise ValueError("DataFrame does not contain 'ref' and 'alt' columns required to filter by VariantInterval. Convert provided interval to GenomicInterval (var_interval.gi()) if you intend to filter only by position.")
    return df_slice


def _sanitize_df(df: pd.DataFrame, interval: GenomicInterval):
    df = df.rename(
        columns={'#chr': 'chrom'}
    )
    if interval is not None:
        df = filter_df_to_interval(df, interval)
    return df


def _parse_interval_type(interval_type: str):
    class_name_dict = {
        'genomic_interval': GenomicInterval,
        'variant_interval': VariantInterval,
    }
    assert interval_type in class_name_dict, "interval_type must be either 'genomic_interval' or 'variant_interval'"
    return class_name_dict[interval_type]


def df_to_intervals(df: pd.DataFrame, interval: GenomicInterval = None, extra_columns=(), verbose=False, interval_type='genomic_interval'):
    interval_type_cls = _parse_interval_type(interval_type)
    df = _sanitize_df(df, interval)
    iterator = df.itertuples(index=False, name=None)
    if verbose:
        iterator = tqdm(iterator, total=len(df), desc="Converting df to intervals")
    get_col_idx = {col: i for i, col in enumerate(df.columns)}
    result = [
        interval_type_cls(
            **{
                col: row[get_col_idx[col]] for col in 
                    interval_type_cls.required_fields + list(extra_columns)
            }
        )
        for row in iterator
    ]
    return result


def df_to_genomic_intervals(df: pd.DataFrame, interval: GenomicInterval=None, extra_columns=()):
    return df_to_intervals(
        df,
        interval=interval,
        extra_columns=extra_columns,
        interval_type='genomic_interval'
    )


def df_to_variant_intervals(df: pd.DataFrame, interval: GenomicInterval=None, extra_columns=()):
    return df_to_intervals(
        df,
        interval=interval,
        extra_columns=extra_columns,
        interval_type='variant_interval'
    )

def intervals_to_df(intervals: List[GenomicInterval]):
    assert len(intervals) > 0, "No intervals to convert"
    interval_type_cls = type(intervals[0])
    
    req_fields = interval_type_cls.required_fields
    extra_fields = intervals[0].extra_fields
    fields = req_fields + extra_fields

    result = []
    for interval in intervals:
        assert isinstance(interval, interval_type_cls), "Not all intervals are of the same type!"

        assert interval.extra_fields == extra_fields, "Not all intervals have the same extra fields!"

        result.append(
            [getattr(interval, field) for field in fields]
        )
    return pd.DataFrame(
        result,
        columns=fields
    )


def genomic_intervals_to_df(genomic_intervals: List[GenomicInterval]):
    """
    Defunc, exists as legacy alias of interval_to_df
    """
    return intervals_to_df(genomic_intervals)
