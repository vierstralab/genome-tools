import numpy as np
import pandas as pd

from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from genome_tools.plotting.utils import format_axes_to_interval

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.utils import DataBundle

from genome_tools.plotting.modular_plot.loaders.footprint import (
    PosteriorLoader,
    ProtectedNucleotidesLoader,
    FootprintsIndexLoader,
    DifferentialFootprintLoader,
    FootprintsDataLoader,
    SequenceWeightsFromProtectedNucleotidesLoader,
)
from genome_tools.plotting.modular_plot.loaders.basic import AnnotationRegionsLoader
from genome_tools.plotting.modular_plot.loaders.variant import VariantGenotypeLoader, GroupsByGenotypeLoader
from genome_tools.plotting.modular_plot.loaders.sequence import (
    FastaLoader,
    MotifHitsLoader,
    MotifHitsSelectorLoader
)
from genome_tools.plotting.modular_plot.plot_components.sequence import MotifHitsComponent, SequencePlotComponent

from .basic import SegmentPlotComponent


@uses_loaders(FootprintsIndexLoader)
class FootprintsIndexComponent(SegmentPlotComponent):
    __intervals_attr__ = FootprintsIndexLoader.__intervals_attr__

    def _plot(self, data, ax, **kwargs):
        super()._plot(data, ax, **kwargs)
        ax.set_ylabel('Footprint\nindex')
        return ax

@uses_loaders(PosteriorLoader)
class PosteriorHeatmapComponent(IntervalPlotComponent):

    def _plot(self, data, ax, hspace=0.05, **kwargs):
        """
        main plot function of the component
        always accepts data, ax, **kwargs
        kwargs override any fields in init
        """
        grouped_data = pd.DataFrame(data.interval_posterior).groupby(
            data.grouping_column,
            observed=True
        )

        group_names = pd.unique(data.grouping_column)

        grouped_posteriors = [
            (group, grouped_data.get_group(group)) for group in group_names
        ]

        num_groups = len(grouped_posteriors)
        gs = gridspec.GridSpecFromSubplotSpec(
            num_groups, 1,
            height_ratios=[df.shape[0] for _, df in grouped_posteriors],
            subplot_spec=ax.get_subplotspec(),
            hspace=hspace
        )
        
        row = 0
        for group, df in grouped_posteriors:
            fig = ax.get_figure()
            ax1 = fig.add_subplot(gs[row, :])
            ax1.pcolormesh(1 - np.exp(-df), cmap="Blues", **kwargs)
            ax1.xaxis.set_visible(False)
            for s in ['top', 'right']:
                ax1.spines[s].set_visible(True)
            
            # Formatting
            ax1.set_ylabel(f"{group} (n={df.shape[0]})", rotation=0, ha="right", va="center", fontsize="medium")
            ax1.set_xticks([])  # Hide x-axis ticks
            ax1.set_yticks([])  # Hide y-axis ticks
            
            row += 1
        ax.axis("off")
        ax.set_xlim(data.interval.start, data.interval.end)

        return ax


TFProtectedNucleotidesComponent = SequencePlotComponent.with_loaders(
    PosteriorLoader, FastaLoader, ProtectedNucleotidesLoader,
    new_class_name="TFProtectedNucleotidesComponent",
)


@uses_loaders(FootprintsDataLoader)
class FootprintTrackComponent(IntervalPlotComponent):
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, smpl_idx=0, color='k', exp_color='C1', lw=0.5, kind='pp', **kwargs):
        xs = self.squarify_array(np.arange(data.obs.shape[1] + 1) + data.interval.start)
        if kind == 'pp':
            try:
                getattr(data, 'pp')
            except AttributeError:
                raise AttributeError("FootprintTrackComponent with kind='pp' requires calc_posteriors=True")

            ax.plot(xs, np.repeat(data.pp[smpl_idx, :], 2), color=color, lw=lw, **kwargs)
        elif kind == 'obs/exp':
            ax.plot(xs, np.repeat(data.obs[smpl_idx, :], 2), color=exp_color, lw=lw, **kwargs)
            ax.plot(xs, np.repeat(data.exp[smpl_idx, :], 2), color=color, lw=lw, **kwargs)
        format_axes_to_interval(ax, data.interval)
        return ax
    
    @staticmethod
    def squarify_array(y):
        return np.concatenate([y[:1], np.repeat(y[1:-1], 2), y[-1:]])


@uses_loaders(VariantGenotypeLoader, GroupsByGenotypeLoader, FootprintsDataLoader, DifferentialFootprintLoader)
class DifferentialFootprintsComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data: DataBundle, ax: plt.Axes, cmap='Spectral', vmin=-5, vmax=5, **kwargs):
        """
        main plot function of the component
        always accepts data, ax, **kwargs
        kwargs override any fields in init
        """
        foldchange = data.lfc
        neglog_pval = data.neglog10_pval

        neglog_pval = np.where(foldchange < 0, -neglog_pval, neglog_pval)
        
        norm = colors.Normalize(vmin=vmin, vmax=vmax)
        mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)

        color = mappable.to_rgba(neglog_pval)
        
        x = np.arange(data.interval.start, data.interval.end) + 0.5 # center bars
        ax.bar(x, foldchange, width=1, color=color)
        ax.axhline(0, ls='dashed', color='k')

        ax.xaxis.set_visible(False)
        max_ylim = np.max(np.abs(foldchange)) * 1.05
        ax.set_ylim(-max_ylim, max_ylim)
        
        cbar = inset_axes(ax, width="2.5%", height="60%", loc='upper left',
                         bbox_to_anchor=(1.05 , -0.2, 1, 1),
                         bbox_transform=ax.transAxes, borderpad=0,)

        plt.colorbar(mappable, cax=cbar, orientation='vertical', label='-log10 p-value')
        
        return ax


FPWeightedMotifHitsComponent = MotifHitsComponent.with_loaders(
    *TFProtectedNucleotidesComponent.__required_loaders__,
    SequenceWeightsFromProtectedNucleotidesLoader,
    *MotifHitsComponent.__required_loaders__,
    new_class_name="FPWeightedMotifHitsComponent",
)
