from genome_tools.plotting.utils import clear_spines
from genome_tools.plotting import signal_plot
from genome_tools.plotting.gene_annotation import gene_annotation_plot
from genome_tools.plotting.ideogram import ideogram_plot
from genome_tools.plotting.utils import clear_spines


from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders

from genome_tools.plotting.modular_plot.loaders.basic import IdeogramLoader, GencodeLoader, SignalLoader, SegmentsLoader, SegmentsTabixLoader


@uses_loaders(IdeogramLoader)
class IdeogramComponent(IntervalPlotComponent):
    """Render a chromosome ideogram at the interval position.

    Loaders: `IdeogramLoader`

    Required loader args:
    - `ideogram_data`: pre-loaded ideogram data (see `genome_tools.plotting.ideogram.read_ideogram`)

    Plot kwargs passed to `genome_tools.plotting.ideogram.ideogram_plot`

    Returns: `matplotlib.axes.Axes`
    """
    
    def _plot(self, data, ax, **kwargs):
        ideogram_plot(data.ideogram_data, data.interval.chrom, pos=data.interval.start, ax=ax, **kwargs)
        return ax


@uses_loaders(GencodeLoader)
class GencodeComponent(IntervalPlotComponent):
    """Plot GENCODE gene annotations overlapping the interval.

    Loaders: `GencodeLoader`

    Required loader args:
    - `gencode_annotation_file`: path to a GENCODE GTF

    Plot kwargs:
    - `gene_symbol_exclude_regex`: regex to exclude labels (default provided)
    - additional kwargs are passed to `genome_tools.plotting.gene_annotation.gene_annotation_plot`

    Returns: `matplotlib.axes.Axes`
    """
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
        except ValueError:
            self.logger.warning("No gene annotations found for the interval.")
        ax.set_xticks([])
        ax.set_yticks([])
        clear_spines(ax)
        return ax


@uses_loaders(SignalLoader)
class TrackComponent(IntervalPlotComponent):
    """Plot a continuous signal track over an interval.

    Loaders: `SignalLoader`

    Required loader args (via `SignalLoader`):
    - `signal_file`: path to a BigWig file

    Plot kwargs passed to `genome_tools.plotting.signal_plot`

    Returns: `matplotlib.axes.Axes`
    """
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        signal_plot(data.interval, data.signal, ax=ax, **kwargs)
        return ax


@uses_loaders(SegmentsLoader)
class SegmentBedComponent(IntervalPlotComponent):
    """Plot a continuous signal track over an interval.

    Loaders: `SegmentsLoader`
    Required loader args (via `SegmentsLoader`):
    - `signal_file`: path to a BigWig file

    Plot kwargs passed to `genome_tools.plotting.signal_plot`

    Returns: `matplotlib.axes.Axes`
    """
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        ax.set_xlim(data.interval.start, data.interval.end)
        signal_plot(data.interval, data.signal, ax=ax, **kwargs)
        return ax

SegmentsTabixComponent = SegmentBedComponent.with_loaders(SegmentsTabixLoader, new_class_name='SegmentsTabixComponent')
