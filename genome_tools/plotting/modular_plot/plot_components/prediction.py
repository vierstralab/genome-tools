from genome_tools.plotting.modular_plot.loaders.prediction import AttributionsLoader, IntervalDatasetLoader, DHSDatasetLoader, PredictedSignalLoader, BatchLoader

from genome_tools.plotting.modular_plot.plot_components.sequence import SequencePlotComponent, MotifHitsComponent

from genome_tools.plotting.modular_plot.plot_components.basic import TrackComponent
from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders


# TODO fix other components
@uses_loaders(DHSDatasetLoader, AttributionsLoader)
class AttributionsComponent(SequencePlotComponent):
    
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        ax = super()._plot(data, ax, **kwargs)
        ax.axhline(0, color='black', lw=0.25, ls='--')
        return ax



AttributionsComponentFromBatch = AttributionsComponent.with_loaders(
    BatchLoader, AttributionsLoader,
    new_class_name='AttributionsComponentFromBatch',
)


AttributionsWeightedMotifHitsComponent = MotifHitsComponent.with_loaders(
    DHSDatasetLoader, AttributionsLoader, *MotifHitsComponent.__required_loaders__,
    new_class_name='AttributionsWeightedMotifHitsComponent',
)


PredictedSignalComponent = TrackComponent.with_loaders(
    IntervalDatasetLoader, PredictedSignalLoader,
    new_class_name='PredictedSignalComponent',
)
