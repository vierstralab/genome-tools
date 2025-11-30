from genome_tools.plotting.modular_plot.loaders.prediction import AttributionsLoader

from genome_tools.plotting.modular_plot.plot_components.sequence import SequencePlotComponent
from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders


@uses_loaders(AttributionsLoader)
class AttributionsComponent(SequencePlotComponent):
    
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        ax = super()._plot(data, ax, **kwargs)
        ax.axhline(0, color='black', lw=0.25, ls='--')
        return ax