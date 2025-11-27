from typing import List
from genome_tools import GenomicInterval


from genome_tools.plotting.sequence import plot_motif_logo

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.loaders.sequence import MotifHitsLoader, MotifHitsSelectorLoader
from genome_tools.plotting.modular_plot.loaders.basic import AnnotationRegionsLoader


# TODO: add sequence plot

@uses_loaders(AnnotationRegionsLoader, MotifHitsLoader, MotifHitsSelectorLoader)
class MotifHitsComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        ax.axis('off')
        axes = self.add_axes_at_intervals(data.motif_intervals, data.interval, ax=ax)
        self.plot_motifs_for_intervals(data.motif_intervals, axes)
        return ax, axes
    
    @staticmethod
    def plot_motifs_for_intervals(motif_intervals: List[GenomicInterval], axes):
        assert len(motif_intervals) == len(axes)
        for interval, ax in zip(motif_intervals, axes):
            plot_motif_logo(interval.pfm_matrix, rc=interval.orient == '-', font='IBM Plex Mono', ax=ax)
            ax.set_xlabel(interval.tf_name, labelpad=0.5)
        return axes

