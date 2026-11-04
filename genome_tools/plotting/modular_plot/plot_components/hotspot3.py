import numpy as np
import pandas as pd

import matplotlib.pyplot as plt


from genome_tools.plotting.utils import format_axes_to_interval
from genome_tools.plotting.modular_plot import IntervalPlotComponent, uses_loaders
from genome_tools.plotting.modular_plot.utils import DataBundle
from genome_tools.plotting import signal_plot

from genome_tools.plotting.modular_plot.loaders.hotspot3 import AggCutcountsLoader, PerBpBackgroundTrackLoader, HighSignalMaskLoader, SmoothedSignalLoader, BackgroundDensityLoader, BackgroundBigwigLoader
from genome_tools.plotting.modular_plot.loaders.basic import ParquetSignalLoader, SignalLoader

from hotspot3.peak_calling import find_stretches

from .basic import TrackComponent


@uses_loaders(AggCutcountsLoader, PerBpBackgroundTrackLoader, HighSignalMaskLoader)
class SignalAndMeanBGComponent(IntervalPlotComponent):
    """Per-base cut counts aggregated as in hotspot3, colored where above the fitted background
    threshold (`hs_color`) and gray elsewhere (`bg_color`), with the threshold as a dotted, shaded line.

    Loaders: `AggCutcountsLoader`, `PerBpBackgroundTrackLoader`, `HighSignalMaskLoader`

    Required loader args:
    - `cutcounts_tabix_file`: `{id}.cutcounts.bed.gz`
    - `chrom_sizes`: `pd.DataFrame` indexed by chromosome with a `len` column
    - `parquet_bg_file`: `debug/{id}.fit_params.parquet`
    Optional loader args:
    - `mappable_bases_tabix`: mappable bases file passed to hotspot3 (default None, all bases mappable)

    Plot kwargs:
    - `stride`, `fit_threshold_stride`: plot every n-th base of the signal / of the threshold
      (default 20). Larger values plot large regions (e.g. whole chromosomes) faster; use 1 to
      plot every base of a short region.
    - `hs_color`, `bg_color`, `line_lw`, `fit_threshold_alpha`

    Returns: `matplotlib.axes.Axes`
    """
    def _plot(self, data: DataBundle, ax: plt.Axes, stride=20, hs_color='#31a354', bg_color='#C0C0C0', line_lw=0.25, fit_threshold_stride=20, fit_threshold_alpha=0.1, **kwargs):
        
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
        fit_threshold = data.fit_threshold[::fit_threshold_stride]
        x = np.linspace(*xlim, len(fit_threshold))
        ax.plot(x, fit_threshold, color='grey', ls='dotted', lw=line_lw)
        ax.fill_between(
            x,
            np.zeros_like(fit_threshold),
            fit_threshold,
            alpha=fit_threshold_alpha,
            color='#3f5299',
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


DensityAndBackgroundComponent = SignalAndMeanBGComponent.with_loaders(
    SignalLoader,
    BackgroundBigwigLoader,
    HighSignalMaskLoader,
    new_class_name='DensityAndBackgroundComponent',
    doc="""Same plot as `SignalAndMeanBGComponent` from hotspot3 bigWigs instead of cut counts and
    debug parquet files: normalized density, colored above the background, with the background
    as the dotted line.

    Loaders: `SignalLoader`, `BackgroundBigwigLoader`, `HighSignalMaskLoader`

    Required loader args:
    - `signal_file`: `{id}.normalized_density.bw`
    - `background_bw`: `{id}.background.bw`

    Plot kwargs as in `SignalAndMeanBGComponent`, e.g. `stride` and `fit_threshold_stride`: plot
    every n-th base of the density / of the background (default 20). Raise them for large regions
    to plot faster; use 1 to plot every base.

    Example:
        DensityAndBackgroundComponent(signal_file='AG95384.normalized_density.bw',
                                      background_bw='AG95384.background.bw')

    Returns: `matplotlib.axes.Axes`
    """
)


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


@uses_loaders(SignalLoader, BackgroundDensityLoader)
class SignalNoBackgroundComponent(IntervalPlotComponent):
    def _plot(self, data, ax, **kwargs):
        plot_data = np.clip(data.signal - data.background_density, a_min=0, a_max=None)
        signal_plot(data.interval, plot_data, ax=ax, **kwargs)
        return ax
