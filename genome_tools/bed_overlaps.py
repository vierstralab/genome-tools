import pandas as pd
import pyranges as pr


def distance_to_closest(ref_df: pd.DataFrame, map_df: pd.DataFrame):
    """
    For each row of ref_df, distance (bp) to the closest interval in map_df.
    0 if overlapping, NaN if map_df has no intervals on that chromosome.
    Distances follow pyranges/bedtools convention: book-ended intervals are 1 bp apart.

    Requires pyranges (pip install "genome-tools[pyranges]").

    Returns:
        pd.Series aligned to ref_df.index
    """

    coords = {'#chr': 'Chromosome', 'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}

    def to_pyranges(df):
        df = df[[c for c in coords if c in df.columns]].rename(columns=coords)
        assert list(df.columns) == ['Chromosome', 'Start', 'End'], \
            "DataFrame must have one of '#chr'/'chrom' and 'start', 'end' columns"
        return df.reset_index(drop=True)

    ref = to_pyranges(ref_df)
    ref['_ref_idx'] = ref.index

    nearest = pr.PyRanges(ref).nearest(pr.PyRanges(to_pyranges(map_df))).df
    # ties give several rows per ref interval with the same distance
    dist = nearest.groupby('_ref_idx')['Distance'].min()
    return pd.Series(
        dist.reindex(ref.index).to_numpy(),
        index=ref_df.index,
        name='distance',
    )
