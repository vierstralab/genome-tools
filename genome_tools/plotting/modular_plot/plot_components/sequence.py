import numpy as np
from typing import List
from genome_tools import GenomicInterval


from genome_tools.plotting.utils import add_axes_at_intervals, pack_rows
from genome_tools.plotting.sequence import plot_motif_logo, seq_plot

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.loaders.sequence import MotifHitsLoader, MotifHitsSelectorLoader
from genome_tools.plotting.modular_plot.loaders.sequence import FastaLoader, OHESequenceLoader
from genome_tools.plotting.modular_plot.loaders.basic import AnnotationRegionsLoader

from genome_tools.plotting.sequence import seq_plot


@uses_loaders(FastaLoader, OHESequenceLoader)
class SequencePlotComponent(IntervalPlotComponent):

    def _plot(self, data, ax, **kwargs):
        ax.axis('off')
        seq_plot(data.matrix, ax=ax, offset=data.interval.start, **kwargs)
        return ax


@uses_loaders(AnnotationRegionsLoader, MotifHitsLoader, MotifHitsSelectorLoader)
class MotifHitsComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, pack=False, pad=1):
        ax.axis('off')
        intervals = data.motif_intervals
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
