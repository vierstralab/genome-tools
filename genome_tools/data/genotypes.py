"""Readers for genotyping QC outputs (plink2 KING relatedness, bcftools stats)
and relatedness-based clustering of individuals."""

import re
from io import StringIO

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy


def read_king_matrix(plink_prefix, clip_cutoff=0.4):
    """
    Read a square KING relatedness matrix written by `plink2 --make-king square`.

    Parameters
    ----------
    plink_prefix : str
        Path to the plink2 output without the `.king`/`.king.id` suffix.
    clip_cutoff : float
        Kinship coefficients below this value (and NaNs) are set to 0.

    Returns
    -------
    pd.DataFrame
        Individuals x individuals kinship matrix.
    """
    indivs = np.loadtxt(f'{plink_prefix}.king.id', dtype=str)
    rel_mat = np.loadtxt(f'{plink_prefix}.king')
    mat = pd.DataFrame(rel_mat, index=indivs, columns=indivs)
    mat[mat < clip_cutoff] = 0
    mat[np.isnan(mat)] = 0
    return mat


def read_bcftools_psc(bcftools_stats):
    """
    Read the per-sample counts (PSC) section of `bcftools stats` output.

    Parameters
    ----------
    bcftools_stats : str
        Path to the `bcftools stats` output. A file with only the PSC lines
        (e.g. from `grep PSC`) works as well.

    Returns
    -------
    pd.DataFrame
        One row per sample; the `sample` column is renamed to `indiv_id`.
    """
    header = None
    rows = []
    with open(bcftools_stats) as f:
        for line in f:
            if line.startswith('# PSC\t'):
                header = line
            elif line.startswith('PSC\t'):
                rows.append(line)
    if header is None:
        raise ValueError(f'No PSC header found in {bcftools_stats}')

    stats = pd.read_table(StringIO(header + ''.join(rows)))
    stats = stats.iloc[:, 2:]
    stats.columns = [re.sub(r'^\[\d+\]', '', c) for c in stats.columns]
    return stats.rename(columns={'sample': 'indiv_id'})


def cluster_relatedness(mat, linkage=None, max_dist=0.1, id_prefix='INDIV_'):
    """
    Group samples of the same individual by clustering a kinship matrix.

    Parameters
    ----------
    mat : pd.DataFrame
        Square kinship matrix, e.g. from `read_king_matrix`.
    linkage : np.ndarray, optional
        Linkage matrix of `mat` with rows and columns sorted by sample id. If not
        provided, complete-linkage clustering of the rows by correlation.
    max_dist : float
        Distance threshold for `scipy.cluster.hierarchy.fcluster`.
    id_prefix : str
        Prefix of the new individual ids. Clusters are numbered by size (largest
        first), ties broken by the smallest sample id in the cluster, so the ids
        don't depend on the order of samples in `mat`.

    Returns
    -------
    clusters : pd.DataFrame
        Columns `indiv_id` (original id), `cluster_id`, `genotype_cluster`,
        `count` (cluster size) and `clustered_indiv_id` (new individual id),
        one row per sample in the order of `mat`.
    linkage : np.ndarray
        Linkage matrix of the samples sorted by id, can be passed to
        `genome_tools.plotting.relatedness.relatedness_heatmap`.
    """
    # cluster samples sorted by id, so the result doesn't depend on the order in mat
    sorted_ids = mat.index.sort_values()
    sorted_mat = mat.loc[sorted_ids, sorted_ids]
    if linkage is None:
        linkage = hierarchy.linkage(sorted_mat, method='complete', metric='correlation')
    cl = hierarchy.fcluster(linkage, max_dist, criterion='distance')

    clusters = pd.DataFrame({'indiv_id': sorted_mat.index, 'cluster_id': cl})

    sizes = clusters.groupby('cluster_id').agg(
        count=('indiv_id', 'size'),
        min_indiv_id=('indiv_id', 'min'),
    ).sort_values(['count', 'min_indiv_id'], ascending=[False, True])
    sizes['genotype_cluster'] = np.arange(len(sizes))

    clusters = clusters.merge(
        sizes[['genotype_cluster', 'count']].reset_index(), on='cluster_id'
    )
    clusters['clustered_indiv_id'] = id_prefix + (clusters['genotype_cluster'] + 1).astype(str).str.zfill(4)

    clusters = pd.DataFrame({'indiv_id': mat.index}).merge(clusters, on='indiv_id', how='left')
    return clusters, linkage
