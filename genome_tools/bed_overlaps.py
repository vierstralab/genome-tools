import pandas as pd
import pyranges as pr
import numpy as np
from tqdm.auto import tqdm


def to_pyranges(df: pd.DataFrame):

    coords = {'#chr': 'Chromosome', 'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}

    df = df[[c for c in coords if c in df.columns]].rename(columns=coords)
    assert list(df.columns) == ['Chromosome', 'Start', 'End'], \
        "DataFrame must have one of '#chr'/'chrom' and 'start', 'end' columns"
    return df.reset_index(drop=True)


def _overlap_indicator(ref_pr: pr.PyRanges, map_df: pd.DataFrame, n: int) -> np.ndarray:
    hits = ref_pr.overlap(pr.PyRanges(to_pyranges(map_df))).df
    out = np.zeros(n, dtype=bool)
    if '_ref_idx' in hits.columns:  # no overlaps gives a DataFrame without columns
        out[hits['_ref_idx'].to_numpy()] = True
    return out


def overlaps_any(ref_df: pd.DataFrame, map_df: pd.DataFrame = None, **map_dfs: pd.DataFrame):
    """
    For each row of ref_df, True if it overlaps at least one interval in map_df(s), else False.
    Uses bed coordinates (half-open).

    Usage:
        overlaps_any(peaks, tss_df)                  -> pd.Series (bool, name='overlap')
        overlaps_any(peaks, tss=tss_df, enh=enh_df)  -> pd.DataFrame, one bool column per kwarg

    Pass either map_df or keyword maps, not both.
    Raises ValueError if no map is provided.

    Requires pyranges (pip install "genome-tools[pyranges]").
    """
    if map_df is None and not map_dfs:
        raise ValueError(
            "No map DataFrame(s) provided. Pass either a single positional map_df "
            "or one or more keyword map DataFrames (e.g. tss=tss_df, enh=enh_df)."
        )
    if map_df is not None and map_dfs:
        raise ValueError(
            "Both positional map_df and keyword map DataFrames "
            f"({', '.join(map_dfs)}) were provided. Please use only one of them."
        )
    ref = to_pyranges(ref_df)
    ref['_ref_idx'] = ref.index
    ref_pr = pr.PyRanges(ref)  # built once, reused for every map
    n = len(ref)

    if map_dfs:
        out = {}
        for name, m in tqdm(map_dfs.items(), total=len(map_dfs), desc='overlaps_any'):
            out[name] = _overlap_indicator(ref_pr, m, n)
        return pd.DataFrame(out, index=ref_df.index)

    return pd.Series(_overlap_indicator(ref_pr, map_df, n), index=ref_df.index, name='overlap')


def distance_to_closest(ref_df: pd.DataFrame, map_df: pd.DataFrame):
    """
    For each row of ref_df, distance (bp) to the closest interval in map_df.
    0 if overlapping, NaN if map_df has no intervals on that chromosome.
    Distances follow pyranges/bedtools convention: book-ended intervals are 1 bp apart.

    Requires pyranges (pip install "genome-tools[pyranges]").

    Returns:
        pd.Series aligned to ref_df.index
    """
    ref = to_pyranges(ref_df)
    ref['_ref_idx'] = ref.index

    nearest = pr.PyRanges(ref).nearest(pr.PyRanges(to_pyranges(map_df))).df
    if '_ref_idx' not in nearest.columns:  # map has no intervals on any ref chromosome
        dist = pd.Series(dtype=float)
    else:
        # ties give several rows per ref interval with the same distance
        dist = nearest.groupby('_ref_idx')['Distance'].min()
    return pd.Series(
        dist.reindex(ref.index).to_numpy(),
        index=ref_df.index,
        name='distance',
    )
