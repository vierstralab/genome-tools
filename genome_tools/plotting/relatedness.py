import sys

import matplotlib.pyplot as plt
from scipy.cluster import hierarchy


def relatedness_heatmap(mat, linkage=None, ax=None, ax_dendrogram=None, cmap='Blues', **kwargs):
    """
    Plot a kinship matrix ordered by hierarchical clustering.

    Parameters
    ----------
    mat : pd.DataFrame
        Square kinship matrix, e.g. from `genome_tools.data.genotypes.read_king_matrix`.
    linkage : np.ndarray, optional
        Linkage matrix of the samples sorted by id, e.g. from
        `genome_tools.data.genotypes.cluster_relatedness`. If not provided,
        complete-linkage clustering of the rows by correlation.
    ax : matplotlib.axes.Axes, optional
        Axes for the heatmap.
    ax_dendrogram : matplotlib.axes.Axes, optional
        Axes for the dendrogram; not drawn if not provided.
    **kwargs
        Passed to `ax.imshow`.

    Returns
    -------
    matplotlib.axes.Axes
        The heatmap axes.
    """
    if ax is None:
        ax = plt.gca()
    # linkage from cluster_relatedness refers to samples sorted by id
    sorted_ids = mat.index.sort_values()
    sorted_mat = mat.loc[sorted_ids, sorted_ids]
    if linkage is None:
        linkage = hierarchy.linkage(sorted_mat, method='complete', metric='correlation')

    # dendrogram recursion depth grows with the number of samples
    sys.setrecursionlimit(max(sys.getrecursionlimit(), 10000))
    dendro = hierarchy.dendrogram(
        linkage,
        no_plot=ax_dendrogram is None,
        ax=ax_dendrogram,
        no_labels=True,
    )
    order = dendro['leaves']

    ax.imshow(
        sorted_mat.iloc[order, order].to_numpy(),
        cmap=cmap,
        interpolation='none',
        aspect='equal',
        **kwargs,
    )
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color('k')

    return ax
