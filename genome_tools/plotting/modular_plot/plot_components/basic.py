from genome_tools.plotting import segment_plot
from genome_tools.plotting.utils import clear_spines
from genome_tools.plotting import signal_plot, segment_plot
from genome_tools.plotting.gene_annotation import gene_annotation_plot
from genome_tools.plotting.ideogram import ideogram_plot
from genome_tools.plotting.utils import clear_spines


from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders

from genome_tools.plotting.modular_plot.loaders.basic import IdeogramLoader, GencodeLoader, SignalLoader

class SegmentPlotComponent(IntervalPlotComponent):
    __intervals_attr__ = 'intervals'

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        segment_plot(data.interval, getattr(data, self.__intervals_attr__), ax=ax, **kwargs)
        ax.set_xticks([])
        ax.set_yticks([])
        clear_spines(ax)
        return ax
    

@uses_loaders(SignalLoader)
class TrackComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        ax.set_xlim(data.interval.start, data.interval.end)
        signal_plot(data.interval, data.signal, ax=ax, **kwargs)
        return ax


@uses_loaders(IdeogramLoader)
class IdeogramComponent(IntervalPlotComponent):
    
    def _plot(self, data, ax, **kwargs):
        ideogram_plot(data.ideogram_data, data.interval.chrom, pos=data.interval.start, ax=ax, **kwargs)
        return ax


@uses_loaders(GencodeLoader)
class GencodeComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, gene_symbol_exclude_regex=r'^ENSG|^MIR|^LINC|.*-AS.*', **kwargs):
        try:
            gene_annotation_plot(
                data.interval,
                data.gencode_annotation_file,
                ax=ax,
                gene_symbol_exclude_regex=gene_symbol_exclude_regex,
                **kwargs
            )
            ax.set_xticks([])
            ax.set_yticks([])
            clear_spines(ax)
        except ValueError:
            self.logger.warning("No gene annotations found for the interval.")
        return ax