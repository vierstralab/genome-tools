# This file defines abstract base classes for plot components in a modular plotting framework.
# They can't be used directly but provide common functionality for derived classes.

from genome_tools.plotting.utils import clear_spines
from genome_tools.plotting import segment_plot
from genome_tools.plotting.utils import clear_spines

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from genome_tools.plotting.modular_plot import IntervalPlotComponent


class SegmentPlotComponent(IntervalPlotComponent):
    __intervals_attr__ = 'intervals'
    """Abstract segment renderer for interval-based components.

    Expects the loader to attach a list of `GenomicInterval` objects
    on the `DataBundle` under the attribute named by `__intervals_attr__`.

    Extend this class for specific segment sources (e.g., DHS, peaks).
    Override `__intervals_attr__` when the loader uses a custom attribute.
    Plot kwargs: forwarded to `segment_plot`.

    Returns: `matplotlib.axes.Axes`
    """
    @IntervalPlotComponent.set_xlim_interval
    def _plot(self, data, ax, **kwargs):
        segment_plot(data.interval, getattr(data, self.__intervals_attr__), ax=ax, **kwargs)
        ax.set_xticks([])
        ax.set_yticks([])
        clear_spines(ax)
        return ax
    

class HeatmapPlotComponent(IntervalPlotComponent):

    def _plot(self, data, ax: plt.Axes, hspace=0.05, cmap="Blues", **kwargs):
        """
        main plot function of the component
        always accepts data, ax, **kwargs
        kwargs override any fields in init
        """
        heatmap_data: pd.DataFrame = data.heatmap_data
        grouping_column: pd.Series = data.grouping_column
        grouped_data = heatmap_data.groupby(
            grouping_column,
            observed=True
        )

        group_names = pd.unique(grouping_column)

        grouped_data = [
            (group, grouped_data.get_group(group)) for group in group_names
        ]

        num_groups = len(grouped_data)
        gs = gridspec.GridSpecFromSubplotSpec(
            num_groups, 1,
            height_ratios=[df.shape[0] for _, df in grouped_data],
            subplot_spec=ax.get_subplotspec(),
            hspace=hspace
        )
        axes = []
        row = 0
        for group, df in grouped_data:
            fig = ax.get_figure()
            ax1 = fig.add_subplot(gs[row, :])
            ax1.pcolormesh(df, cmap=cmap, **kwargs)
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