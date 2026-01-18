from typing import List
from genome_tools import GenomicInterval


from genome_tools.plotting.utils import add_axes_at_intervals, pack_rows
from genome_tools.plotting.sequence import plot_motif_logo, seq_plot

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.loaders.sequence import (
    MotifHitsLoader,
    MotifHitsSelectorLoader,
    FastaLoader,
    OHESequenceLoader
)
from genome_tools.plotting.modular_plot.loaders.basic import AnnotationRegionsLoader


# TODO fix other components
@uses_loaders(FastaLoader, OHESequenceLoader)
class SequencePlotComponent(IntervalPlotComponent):
    """Render DNA sequence (one-hot encoded) for the interval.

    Loaders: ``FastaLoader``, ``OHESequenceLoader``

    Required loader args:
    - ``fasta_file``: path to the genome FASTA

    Plot kwargs: forwarded to ``seq_plot`` (e.g., color mapping).

    Returns: ``matplotlib.axes.Axes``
    """

    def _plot(self, data, ax, preserve_aspect_ratio=False, **kwargs):
        ax.axis('off')
        seq_plot(data.matrix, ax=ax, offset=data.interval.start, preserve_aspect_ratio=preserve_aspect_ratio, **kwargs)
        return ax


@uses_loaders(AnnotationRegionsLoader, MotifHitsLoader, MotifHitsSelectorLoader)
class MotifHitsComponent(IntervalPlotComponent):
    """Plot motif logos at selected hits within annotation regions.

    Loaders: ``AnnotationRegionsLoader``, ``MotifHitsLoader``, ``MotifHitsSelectorLoader``

    Required loader args:
    - ``annotation_regions``: list[GenomicInterval]
    - ``motif_annotations_path``: tabix-indexed motif hits file
    - ``motif_meta``: DataFrame indexed by motif_id with metadata

    Selector loader args (examples):
    - ``choose_by``: one of {'dg','weighted_dg','ddg','concordant_ddg'}
    - ``n_top_hits``: int
    - ``motif_hits_threshold``: float | None

    Plot kwargs:
    - ``pack``: bool to pack overlapping rows
    - ``pad``: int padding between packed rows

    Returns: tuple (top axes, list of per-hit axes)
    """

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, pack=False, pad=1):
        ax.axis('off')
        intervals = data.motif_intervals
        if len(intervals) == 0:
            self.logger.warning(f"No motif hits passing selection criteria for intervals {data.annotation_regions}")
            return ax, []
        if pack:
           _, intervals = pack_rows(intervals, pad=pad)
        axes = add_axes_at_intervals(intervals, data.interval, ax=ax)
        self.plot_motifs_for_intervals(intervals, axes)
        return ax, axes
    
    @staticmethod
    def plot_motifs_for_intervals(motif_intervals: List[GenomicInterval], axes):
        assert len(motif_intervals) == len(axes)
        for interval, ax in zip(motif_intervals, axes):
            plot_motif_logo(interval.pfm_matrix, rc=interval.orient == '-', font='IBM Plex Mono', ax=ax)
            ax.set_xlabel(interval.tf_name, labelpad=0.5)
        return axes

