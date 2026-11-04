from genome_tools.plotting.modular_plot.loaders.prediction import (
    BatchLoader,
    BatchFromAnndataLoader,
    BatchFromSteppedIntervalLoader,
    BatchFromIntervalCenterLoader,

    AttributionsLoader,
    AlignedAttributionsLoader, 
    BetweenSpeciesAlignedAttributionsLoader,

    PredictedSignalLoader
)

from genome_tools.plotting.modular_plot.plot_components.sequence import SequencePlotComponent, MotifHitsComponent

from genome_tools.plotting.modular_plot.plot_components.basic import TrackComponent
from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders


# TODO fix other components
PredictedSignalComponent = TrackComponent.with_loaders(
    BatchFromSteppedIntervalLoader, PredictedSignalLoader,
    new_class_name='PredictedSignalComponent',
    doc="""Plot a vinson model's predicted signal (e.g. DNase density) for one sample across the interval.

    The model is run on sequences centered every `step` bp across the interval and
    the predictions are linearly interpolated to per-base resolution. Plotted like
    `TrackComponent`, so it is easy to stack next to the observed track.

    Loaders: `BatchFromSteppedIntervalLoader`, `PredictedSignalLoader`

    Required loader args:
    - `sample_id`: sample to predict for (row of `embeddings`)
    - `embeddings`: DataFrame of sample embeddings indexed by sample_id
      (e.g. `anndata.obsm['motif_embeddings']`)
    - `fasta_file`: path to the genome FASTA
    - `model_wrapper`: `vinson.postprocessing.interpretation.ModelWrapper`

    Optional loader args:
    - `step` (20): distance in bp between prediction centers
    - `genotype_file` (None): VCF to inject sample genotypes into the sequence
    - `interp1d_kind` ('linear'): passed to `scipy.interpolate.interp1d`

    Plot kwargs passed to `genome_tools.plotting.signal_plot`

    Example:
        PredictedSignalComponent(sample_id='AG95451', interval_key='med',
                                 height=0.5, color='tab:orange')
        # model_wrapper, fasta_file and embeddings are usually passed once to IntervalPlotter

    Returns: `matplotlib.axes.Axes`
    """,
)

PredictedSignalFromBatchComponent = TrackComponent.with_loaders(
    BatchLoader, PredictedSignalLoader,
    new_class_name='PredictedSignalFromBatchComponent',
    doc="""Same as `PredictedSignalComponent`, but uses a pre-built `batch`
    (dict with 'ohe_seq', 'embed', 'summit') instead of building one from the interval.

    Loaders: `BatchLoader`, `PredictedSignalLoader`

    Required loader args: `batch`, `model_wrapper`

    Returns: `matplotlib.axes.Axes`
    """,
)


# Attributions for custom region and sample
AttributionsFromRegionComponent = SequencePlotComponent.with_loaders(
    BatchFromIntervalCenterLoader, AttributionsLoader, AlignedAttributionsLoader,
    new_class_name='AttributionsFromRegionComponent',
    doc="""Plot per-nucleotide attributions (sequence logo scaled by attribution)
    of a vinson model for one sample, for a single model input centered on the interval center.

    Loaders: `BatchFromIntervalCenterLoader`, `AttributionsLoader`, `AlignedAttributionsLoader`

    Required loader args: `sample_id`, `embeddings`, `fasta_file`, `model_wrapper`

    Optional loader args: `genotype_file`, `n_shuffles` (20), `random_state` (42),
    `print_convergence_deltas` (False). See `AttributionsFromBatchComponent`.

    Returns: `matplotlib.axes.Axes`
    """,
)

AttributionsWeightedMotifHitsFromRegionComponent = MotifHitsComponent.with_loaders(
    *AttributionsFromRegionComponent.__required_loaders__,
    *MotifHitsComponent.__required_loaders__,
    new_class_name='AttributionsWeightedMotifHitsFromRegionComponent',
    doc="""Motif logos for the best motif hits in `annotation_regions`, using
    attributions computed as in `AttributionsFromRegionComponent`.
    See `AttributionsWeightedMotifHitsFromBatchComponent` for arguments.

    Returns: tuple (top axes, list of per-hit axes)
    """,
)


# Attributions from custom batch
AttributionsFromBatchComponent = SequencePlotComponent.with_loaders(
    BatchLoader, AttributionsLoader, AlignedAttributionsLoader,
    new_class_name='AttributionsFromBatchComponent',
    doc="""Plot per-nucleotide attributions of a vinson model as a sequence logo
    (each base scaled by its attribution score).

    Attributions are computed with `model_wrapper.get_sequence_attributions`
    (DeepLIFT/SHAP with `n_shuffles` dinucleotide-shuffled references) for a
    single-example `batch` and re-aligned to the component's interval, which must
    overlap the model input. Use a small interval (tens of bp) to keep letters readable.

    Loaders: `BatchLoader`, `AttributionsLoader`, `AlignedAttributionsLoader`

    Required loader args:
    - `batch`: single-example batch dict with 'ohe_seq', 'embed', 'chrom', 'summit'
      (e.g. from `vinson.inference.assemble_batch`)
    - `model_wrapper`: `vinson.postprocessing.interpretation.ModelWrapper`

    Optional loader args:
    - `n_shuffles` (20): number of shuffled references; more is slower but more stable
    - `random_state` (42)
    - `print_convergence_deltas` (False)

    Plot kwargs passed to `genome_tools.plotting.sequence.seq_plot`

    Sets `data.matrix` (len(interval) x 4) and `data.sequence_weights`
    (per-base attribution sums), which motif components reuse.

    Example:
        AttributionsFromBatchComponent(interval_key='small', height=0.06, n_shuffles=100)

    Returns: `matplotlib.axes.Axes`
    """,
)

AttributionsWeightedMotifHitsFromBatchComponent = MotifHitsComponent.with_loaders(
    *AttributionsFromBatchComponent.__required_loaders__,
    *MotifHitsComponent.__required_loaders__,
    new_class_name='AttributionsWeightedMotifHitsFromBatchComponent',
    doc="""Plot motif logos for motif hits inside `annotation_regions`, computed
    from the same attributions as `AttributionsFromBatchComponent`.

    Motif hits are read from a tabix-indexed motif scan; for each annotation region the
    best hit is kept (or all hits above `motif_hits_threshold`). Use
    `choose_by='weighted_dg'` to rank hits by dG weighted with the attributions.

    Loaders: `BatchLoader`, `AttributionsLoader`, `AlignedAttributionsLoader`,
    `AnnotationRegionsLoader`, `MotifHitsLoader`, `MotifHitsSelectorLoader`

    Required loader args:
    - `batch`, `model_wrapper`: as in `AttributionsFromBatchComponent`
    - `annotation_regions`: list of `GenomicInterval` that overlap the interval
    - `motif_annotations_path`: tabix-indexed motif hits (MOODS output) with columns
      '#chr', 'start', 'end', 'motif_id', 'dg', 'orient', 'seq'
    - `motif_meta`: DataFrame indexed by motif_id with motif metadata
      (needs 'pfm' (path to .pfm file) and 'tf_name' columns)

    Optional loader args:
    - `choose_by` ('dg'): one of 'dg', 'weighted_dg', 'ddg', 'concordant_ddg'
    - `n_top_hits` (1), `motif_hits_threshold` (None), `min_motif_overlap` (0.9)
    - attribution args of `AttributionsFromBatchComponent`

    Plot kwargs:
    - `pack` (False): pack overlapping logos into rows; `pad` (1) between rows
    - other kwargs are passed to `genome_tools.plotting.sequence.plot_motif_logo`

    Example:
        AttributionsWeightedMotifHitsFromBatchComponent(
            interval_key='small', height=0.06, n_shuffles=100,
            motif_annotations_path=motif_annotations_path,
        )
        # annotation_regions, batch, motif_meta are passed to IntervalPlotter or get_interval_data

    Returns: tuple (top axes, list of per-hit axes)
    """,
)


# Attributions from custom batch with between species alignment
BetweenSpeciesAlignedAttributionsFromBatchComponent = SequencePlotComponent.with_loaders(
    BatchLoader, AttributionsLoader, BetweenSpeciesAlignedAttributionsLoader,
    new_class_name='BetweenSpeciesAlignedAttributionsFromBatchComponent',
    doc="""Same as `AttributionsFromBatchComponent`, but the attributions are mapped onto
    the interval through a `BetweenSpeciesMap` (`mapping` loader arg), e.g. to show
    attributions computed on a mouse sequence in human coordinates.

    Returns: `matplotlib.axes.Axes`
    """,
)

BetweenSpeciesAlignedAttributionsWeightedMotifHitsFromBatchComponent = MotifHitsComponent.with_loaders(
    *BetweenSpeciesAlignedAttributionsFromBatchComponent.__required_loaders__,
    *MotifHitsComponent.__required_loaders__,
    new_class_name='BetweenSpeciesAlignedAttributionsWeightedMotifHitsFromBatchComponent',
    doc="""Motif hits counterpart of `BetweenSpeciesAlignedAttributionsFromBatchComponent`.
    See `AttributionsWeightedMotifHitsFromBatchComponent` for arguments.

    Returns: tuple (top axes, list of per-hit axes)
    """,
)


# Attributions for DHS from anndata. Defunc
@uses_loaders(BatchFromAnndataLoader, AttributionsLoader, AlignedAttributionsLoader)
class AttributionsFromAnndataComponent(SequencePlotComponent):
    """Defunct. Attributions for a DHS (`dhs_id`) and sample taken from a backed AnnData."""

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        ax = super()._plot(data, ax, **kwargs)
        ax.axhline(0, color='black', lw=0.25, ls='--')
        return ax

AttributionsFromAnndataWeightedMotifHitsComponent = MotifHitsComponent.with_loaders(
    *AttributionsFromAnndataComponent.__required_loaders__, *MotifHitsComponent.__required_loaders__,
    new_class_name='AttributionsFromAnndataWeightedMotifHitsComponent',
    doc="""Defunct. Motif hits counterpart of `AttributionsFromAnndataComponent`.""",
)