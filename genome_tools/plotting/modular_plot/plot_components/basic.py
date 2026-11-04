from genome_tools.plotting.utils import clear_spines
from genome_tools.plotting import signal_plot
from genome_tools.plotting.gene_annotation import gene_annotation_plot
from genome_tools.plotting.ideogram import ideogram_plot
from genome_tools.plotting.utils import clear_spines


from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders

from genome_tools.plotting.modular_plot.loaders.basic import IdeogramLoader, GencodeLoader, SignalLoader, FigwigSignalLoader, SegmentsLoader, SegmentsTabixLoader

from .abstract import SegmentPlotComponent


class EmptyComponent(IntervalPlotComponent):
    """
    Placeholder class to use as Plot component without modifying the ax
    
    Returns: `matplotlib.axes.Axes`
    """
    def _plot(self, data, ax, **kwargs):
        return ax

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
    """Plot a continuous per-base signal track (e.g. DNase density) over an interval.

    Loaders: `SignalLoader`

    Required loader args (via `SignalLoader`):
    - `signal_file`: path to a BigWig file. NaNs are replaced with 0.

    Plot kwargs passed to `genome_tools.plotting.signal_plot`, e.g.
    `color`, `lw`, `fill`. Set y-limits on the returned axes.

    Subclasses that keep this `_plot` and only swap loaders
    (`FigwigTrackComponent`, `PredictedSignalComponent`) need a loader that sets
    `data.signal` to a 1D array of length `len(interval)`.

    Example:
        TrackComponent(signal_file='sample.normalized_density.bw',
                       height=0.5, margins=(0.05, 0.1), color='tab:blue')

    Returns: `matplotlib.axes.Axes`
    """
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        signal_plot(data.interval, data.signal, ax=ax, **kwargs)
        return ax


FigwigTrackComponent = TrackComponent.with_loaders(
    FigwigSignalLoader,
    new_class_name='FigwigTrackComponent',
    doc="""Plot one track of a figwig (.fw) file over an interval.

    Same as `TrackComponent`, but the signal is read from a multi-track figwig
    file (https://github.com/jmschrei/figwig, requires python >= 3.10).

    Loaders: `FigwigSignalLoader`

    Required loader args:
    - `figwig_file`: path to the .fw file
    - `track`: name of the track to plot

    Plot kwargs passed to `genome_tools.plotting.signal_plot`

    Example:
        FigwigTrackComponent(figwig_file='samples.fw', track='AG80080',
                             height=0.5, color='tab:blue')

    Returns: `matplotlib.axes.Axes`
    """
)


SegmentBedComponent = SegmentPlotComponent.with_loaders(SegmentsLoader, new_class_name='SegmentBedComponent')


SegmentsTabixComponent = SegmentPlotComponent.with_loaders(SegmentsTabixLoader, new_class_name='SegmentsTabixComponent')
