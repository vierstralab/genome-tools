import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


from genome_tools.plotting.utils import format_axes_to_interval
from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.utils import DataBundle

from genome_tools.plotting.modular_plot.loaders.hotspot3 import AggCutcountsLoader, PerBpBackgroundTrackLoader, HighSignalMaskLoader, SmoothedSignalLoader
from genome_tools.plotting.modular_plot.loaders.basic import ParquetSignalLoader

from hotspot3.peak_calling import find_stretches

from .basic import TrackComponent


@uses_loaders(AggCutcountsLoader, PerBpBackgroundTrackLoader, HighSignalMaskLoader)
class SignalAndMeanBGComponent(IntervalPlotComponent):
    def _plot(self, data: DataBundle, ax: plt.Axes, stride=500, hs_color='#31a354', bg_color='#C0C0C0', **kwargs):
        
        self.plot_bg_and_signal(
            data.signal,
            data.high_signal_mask,
            data.interval,
            ax=ax,
            stride=stride,
            hs_color=hs_color,
            bg_color=bg_color,
            **kwargs
        )
        xlim = ax.get_xlim()
        ax.plot(np.linspace(*xlim, len(data.fit_threshold)), data.fit_threshold, color='grey', ls='dotted', lw=0.5)
        ax.fill_between(
            np.linspace(*xlim, len(data.fit_threshold)),
            np.zeros_like(data.fit_threshold),
            data.fit_threshold,
            alpha=0.1,
            color='#3f5299',
            rasterized=True,
            lw=0
        )
        
        return ax

    @staticmethod
    def plot_bg_and_signal(data, high_signal_mask, region, ax=None, stride=1, bg_color='#3f5299', hs_color='#C0C0C0', linewidth=0.25, **kwargs):
        if ax is None:
            ax = plt.gca()
        hs_data = np.where(high_signal_mask, data, np.nan)[::stride]
        bg_data = np.where(~high_signal_mask, data, np.nan)[::stride]
        ax.fill_between(np.linspace(region.start, region.end, len(hs_data)), hs_data, linewidth=linewidth, color=hs_color, rasterized=True, **kwargs)
        ax.fill_between(np.linspace(region.start, region.end, len(bg_data)), bg_data, linewidth=linewidth, color=bg_color, rasterized=True, **kwargs)
        format_axes_to_interval(ax, region)
        return ax


@uses_loaders(*SignalAndMeanBGComponent.__required_loaders__, SmoothedSignalLoader)
class SmoothedSignalAndMeanBGComponent(SignalAndMeanBGComponent):
    def _plot(self, data: DataBundle, ax: plt.Axes, **kwargs):
        super()._plot(data, ax=ax, **kwargs)
        xlim = ax.get_xlim()
        ax.plot(
            np.linspace(*xlim, len(data.smoothed_signal)),
            data.smoothed_signal,
            color='#3f5299',
            ls='-',
            lw=0.5
        )
        return ax


@uses_loaders(ParquetSignalLoader)
class FdrComponent(TrackComponent):
    def _plot(self, data: DataBundle, ax: plt.Axes, fdr_tr=0.001, color='k', **kwargs):
        super()._plot(data, ax=ax, color=color, **kwargs)

        log_tr = -np.log10(fdr_tr)
        ax.axhline(log_tr, color=color, lw=0.25, ls='dotted')
        starts, ends = find_stretches(data.signal >= log_tr)
        hl = ax.hlines(
            [log_tr]*len(starts),
            starts + data.interval.start,
            ends + data.interval.start,
            color='#d94f4f',
            lw=0.5,
            zorder=1000,
            antialiased=False
        )
        hl.set_capstyle('round')
        ax.tick_params(which='minor', length=0)
        return ax
