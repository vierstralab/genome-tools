import numpy as np
import pandas as pd

from matplotlib import gridspec
from genome_tools.plotting.utils import format_axes_to_interval

from genome_tools.plotting.sequence import seq_plot

from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders

from genome_tools.plotting.modular_plot.loaders.footprint import PosteriorLoader, FootprintDatasetLoader, ProtectedNucleotidesLoader, FootprintsLoader
from genome_tools.plotting.modular_plot.loaders.sequence import FastaLoader



from .basic import SegmentPlotComponent


@uses_loaders(FootprintsLoader)
class FootprintsComponent(SegmentPlotComponent):
    __intervals_attr__ = FootprintsLoader.__intervals_attr__

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
            subplot_spec=ax,
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


# Just use the usual sequence plot with different loaders attached
@uses_loaders(PosteriorLoader, FastaLoader, ProtectedNucleotidesLoader)
class TFProtectedNucleotidesComponent(IntervalPlotComponent):

    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        """
        main plot function of the component
        always accepts data, ax, **kwargs
        kwargs override any fields in init
        """
        seq_plot(data.matrix, ax=ax, font="IBM Plex Mono", offset=data.interval.start)
        ax.axhline(0, color="black", linewidth=0.1)
        ax.text(-0.02, 0.5, "TF-protected\n nucleotides",
                rotation=0, ha="right", va="center", transform=ax.transAxes, fontsize="medium")
        ax.axis("off")
        return ax



@uses_loaders(FootprintDatasetLoader)
class FootprintTrackComponent(IntervalPlotComponent):
    # DEFUNC FIXME
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, smpl_idx=0, color='k', exp_color='C1', lw=0.5, kind='pp', **kwargs):
        xs = self.squarify_array(np.arange(data.pp.shape[1] + 1) + data.interval.start)
        if kind == 'pp':
            ax.plot(xs, np.repeat(data.pp[smpl_idx, :], 2), color=color, lw=lw, **kwargs)
        elif kind == 'obs/exp':
            ax.plot(xs, np.repeat(data.obs[smpl_idx, :], 2), color=exp_color, lw=lw, **kwargs)
            ax.plot(xs, np.repeat(data.exp[smpl_idx, :], 2), color=color, lw=lw, **kwargs)
        format_axes_to_interval(ax, data.interval)
        return ax
    
    @staticmethod
    def squarify_array(y):
        return np.concatenate([y[:1], np.repeat(y[1:-1], 2), y[-1:]])


