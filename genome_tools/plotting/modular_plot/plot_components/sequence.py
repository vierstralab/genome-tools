from genome_tools.plotting.pwm import plot_motif_logo
from genome_tools.plotting.utils import format_axes_to_interval

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.loaders.sequence import MotifHitsLoader
from genome_tools.plotting.modular_plot.loaders.basic import AnnotationRegionsLoader


@uses_loaders(AnnotationRegionsLoader, MotifHitsLoader)
class MotifComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        ax.axis('off')
        axes = self.add_axes_at_intervals(data.motif_intervals, data.interval, ax=ax)
        self.plot_motifs_for_intervals(data.motif_intervals, axes)
        return ax, axes
    
    @staticmethod
    def plot_motifs_for_intervals(motif_intervals, axes):
        assert len(motif_intervals) == len(axes)
        for interval, ax in zip(motif_intervals, axes):
            plot_motif_logo(interval.pfm, rc=interval.orient == '-', font='IBM Plex Mono', ax=ax)
            format_axes_to_interval(ax, interval)
            ax.set_xlabel(interval.tf_name, labelpad=0.5)

