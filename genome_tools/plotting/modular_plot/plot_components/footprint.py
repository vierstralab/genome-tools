import numpy as np
import pandas as pd

from matplotlib import gridspec
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from genome_tools.plotting.utils import format_axes_to_interval

from genome_tools.plotting.colors.cm import get_vocab_color

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.utils import DataBundle

from genome_tools.plotting.modular_plot.loaders.footprint import (
    PosteriorLoader,
    ProtectedNucleotidesLoader,
    FootprintsIndexLoader,
    DifferentialFootprintLoader,
    FootprintsDataLoader,
    SequenceWeightsFromProtectedNucleotidesLoader,
    GroupsDataLoader
)
from genome_tools.plotting.modular_plot.loaders.variant import VariantGenotypeLoader, GroupsByGenotypeLoader, VariantIntervalLoader

from genome_tools.plotting.modular_plot.plot_components.sequence import FastaLoader, MotifHitsComponent, SequencePlotComponent

from .abstract import SegmentPlotComponent


# TODO fix other components
@uses_loaders(FootprintsIndexLoader)
class FootprintsIndexComponent(SegmentPlotComponent):
    __intervals_attr__ = FootprintsIndexLoader.__intervals_attr__

    def _plot(self, data, ax, **kwargs):
        super()._plot(data, ax, **kwargs)
        ax.set_ylabel('Footprint\nindex')
        return ax


@uses_loaders(PosteriorLoader)
class PosteriorHeatmapComponent(IntervalPlotComponent):

    def _plot(self, data, ax: plt.Axes, hspace=0.05, **kwargs):
        """
        main plot function of the component
        always accepts data, ax, **kwargs
        kwargs override any fields in init
        """
        interval_posterior: pd.DataFrame = data.interval_posterior
        grouping_column: pd.Series = data.grouping_column
        grouped_data = interval_posterior.groupby(
            grouping_column,
            observed=True
        )

        group_names = pd.unique(grouping_column)

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
        axes = []
        row = 0
        for group, df in grouped_posteriors:
            fig = ax.get_figure()
            ax1 = fig.add_subplot(gs[row, :])
            ax1.pcolormesh(1 - np.exp(-df), cmap="Blues", **kwargs)
            ax1.xaxis.set_visible(False)
            for s in ['top', 'right']:
                ax1.spines[s].set_visible(True)
            
            # Formatting
            ax1.set_ylabel(group, rotation=0, ha="right", va="center", fontsize="medium")
            ax1.set_xticks([])  # Hide x-axis ticks
            ax1.set_yticks([])  # Hide y-axis ticks
            axes.append(ax1)
            row += 1
        ax.axis("off")

        return ax, axes


TFProtectedNucleotidesComponent = SequencePlotComponent.with_loaders(
    PosteriorLoader, FastaLoader, ProtectedNucleotidesLoader,
    new_class_name="TFProtectedNucleotidesComponent",
)


@uses_loaders(FootprintsDataLoader)
class FootprintTrackComponent(IntervalPlotComponent):
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data: DataBundle, ax: plt.Axes, sample_ids=None, color='k', lw=0.5, kind='pp', exp_kwargs=None, **kwargs):
        if sample_ids is None:
            sample_ids = data.obs.index.tolist()
        xs = self.squarify_array(np.arange(data.obs.shape[1] + 1) + data.interval.start)
        if kind == 'pp':
            try:
                getattr(data, 'pp')
            except AttributeError:
                raise AttributeError("FootprintTrackComponent with kind='pp' requires calc_posteriors=True")
            plot_data = np.repeat(data.pp.loc[sample_ids, :].mean(axis=0), 2)

        elif kind == 'obs/exp':
            exp_plot_kwargs = kwargs.copy()
            exp_plot_kwargs['color'] = '#C0C0C0'
            exp_plot_kwargs['lw'] = lw
            if exp_kwargs is not None:
                exp_plot_kwargs.update(exp_kwargs)

            ax.fill_between(
                xs,
                np.zeros_like(xs),
                np.repeat(data.exp.loc[sample_ids, :].sum(axis=0), 2),
                **exp_plot_kwargs
            )
            plot_data = np.repeat(data.obs.loc[sample_ids, :].sum(axis=0), 2)
        else:
            raise ValueError(f"Unknown kind '{kind}' for FootprintTrackComponent. Expected 'pp' or 'obs/exp'.")

        ax.fill_between(
            xs,
            np.zeros_like(xs),
            plot_data,
            color=color,
            lw=lw,
            **kwargs
        )
        return ax
    
    @staticmethod
    def squarify_array(y):
        return np.concatenate([y[:1], np.repeat(y[1:-1], 2), y[-1:]])


@uses_loaders(VariantIntervalLoader, VariantGenotypeLoader, GroupsByGenotypeLoader,FootprintsDataLoader)
class AllelicFpTrackComponent(FootprintTrackComponent):

    def _plot(self, data: DataBundle, ax: plt.Axes, kind='obs/exp', vocab='dna', **kwargs):
        groups_data: pd.DataFrame = data.groups_data

        groups = ['AA', 'BB']

        gs = gridspec.GridSpecFromSubplotSpec(
            len(groups),
            1,
            subplot_spec=ax.get_subplotspec(), 
            hspace=0.25
        )
        ax.axis("off")
        fig = ax.get_figure()

        axes = []
        for i, (group, allele) in enumerate(zip(groups, [data.variant_interval.ref, data.variant_interval.alt])):
            color = get_vocab_color(allele, vocab)
            ax_group = fig.add_subplot(gs[i, :])
            format_axes_to_interval(ax_group, data.interval)

            sample_ids = groups_data.query(f'parsed_genotype == "{group}"').index.values
            super()._plot(
                data,
                ax_group,
                sample_ids=sample_ids,
                kind=kind,
                color=color,
                **kwargs
            )
            ax_group.xaxis.set_visible(False)
            axes.append(ax_group)

        return ax, axes


@uses_loaders(GroupsDataLoader, FootprintsDataLoader, DifferentialFootprintLoader)
class DifferentialFootprintsComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data: DataBundle, ax: plt.Axes, cmap='Spectral', vmin=-5, vmax=5):
        """
        main plot function of the component
        always accepts data, ax, **kwargs
        kwargs override any fields in init
        """
        foldchange = data.lfc
        neglog_pval = data.neglog10_pval_windowed

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
        
        cax = ax.inset_axes([1.02, 0.0, 0.025, 0.6], transform=ax.transAxes)

        plt.colorbar(mappable, cax=cax, orientation='vertical', label='-log10 p-value')
        
        return ax, cax


DifferentialFootprintsByGenotypeComponent = DifferentialFootprintsComponent.with_loaders(
    VariantIntervalLoader,
    VariantGenotypeLoader,
    GroupsByGenotypeLoader,
    FootprintsDataLoader, 
    DifferentialFootprintLoader,
    new_class_name='DifferentialFootprintsByGenotypeComponent',
)


FPWeightedMotifHitsComponent = MotifHitsComponent.with_loaders(
    *TFProtectedNucleotidesComponent.__required_loaders__,
    SequenceWeightsFromProtectedNucleotidesLoader,
    *MotifHitsComponent.__required_loaders__,
    new_class_name="FPWeightedMotifHitsComponent",
)
