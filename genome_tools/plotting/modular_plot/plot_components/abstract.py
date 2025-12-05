# This file defines abstract base classes for plot components in a modular plotting framework.
# They can't be used directly but provide common functionality for derived classes.

from genome_tools.plotting.utils import clear_spines
from genome_tools.plotting import segment_plot
from genome_tools.plotting.utils import clear_spines

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
    