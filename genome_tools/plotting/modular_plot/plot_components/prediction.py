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
)


# Attributions for DHS from anndata
@uses_loaders(BatchFromAnndataLoader, AttributionsLoader, AlignedAttributionsLoader)
class AttributionsComponent(SequencePlotComponent):
    
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        ax = super()._plot(data, ax, **kwargs)
        ax.axhline(0, color='black', lw=0.25, ls='--')
        return ax

AttributionsWeightedMotifHitsComponent = MotifHitsComponent.with_loaders(
    *AttributionsComponent.__required_loaders__, *MotifHitsComponent.__required_loaders__,
    new_class_name='AttributionsWeightedMotifHitsComponent',
)


# Attributions for custom region and sample from anndata 
AttributionsFromRegionComponent = AttributionsComponent.with_loaders(
    BatchFromIntervalCenterLoader, AttributionsLoader, AlignedAttributionsLoader,
    new_class_name='AttributionsFromBatchComponent',
)

AttributionsWeightedMotifHitsFromRegionComponent = MotifHitsComponent.with_loaders(
    *AttributionsFromRegionComponent.__required_loaders__,
    *MotifHitsComponent.__required_loaders__,
    new_class_name='AttributionsWeightedMotifHitsFromRegionComponent',
)


# Attributions from custom batch
AttributionsFromBatchComponent = AttributionsComponent.with_loaders(
    BatchLoader, AttributionsLoader, AlignedAttributionsLoader,
    new_class_name='AttributionsFromBatchComponent',
)

AttributionsWeightedMotifHitsFromBatchComponent = MotifHitsComponent.with_loaders(
    *AttributionsFromBatchComponent.__required_loaders__,
    *MotifHitsComponent.__required_loaders__,
    new_class_name='AttributionsWeightedMotifHitsFromBatchComponent',
)


# Attributions from custom batch with between species alignment
BetweenSpeciesAlignedAttributionsFromBatchComponent = AttributionsComponent.with_loaders(
    BatchLoader, AttributionsLoader, BetweenSpeciesAlignedAttributionsLoader,
    new_class_name='BetweenSpeciesAlignedAttributionsFromBatchComponent',
)

BetweenSpeciesAlignedAttributionsWeightedMotifHitsFromBatchComponent = MotifHitsComponent.with_loaders(
    *BetweenSpeciesAlignedAttributionsFromBatchComponent.__required_loaders__,
    *MotifHitsComponent.__required_loaders__,
    new_class_name='BetweenSpeciesAlignedAttributionsWeightedMotifHitsFromBatchComponent',
)
