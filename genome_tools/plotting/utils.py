
import numpy as np

import matplotlib.ticker as mticker

__all__ = ["clear_spines", "get_genomic_axis_tick_formatter", "format_axes_to_interval", "rescale_data"]

def clear_spines(ax):
    """Clear spines from plot"""

    for loc, spine in ax.spines.items():
        spine.set_visible(False)

    ax.patch.set_visible(False)


def get_genomic_axis_tick_formatter(interval):
    """Return a tick formattor for the genomic coordinate axis"""

    if len(interval) > 1e5:
        fmt='{:0.2f} Mb'
        denom=1e6
    else:
        fmt='{:,.0f}'
        denom=1

    def _func(x, pos):
        return fmt.format(x/denom) 
    
    return mticker.FuncFormatter(_func)


def format_axes_to_interval(ax, interval, **kwargs):
    """Format axis appearance"""

    ax.set_xlim(interval.start, interval.end)

    ax.tick_params(
        axis='both',
        direction='out'
    )

    ax.xaxis.set(
        major_locator=mticker.MaxNLocator(3, prune='both'),
        minor_locator=mticker.AutoMinorLocator(4),
        major_formatter=get_genomic_axis_tick_formatter(interval),
    )

    ax.yaxis.set(
        major_locator=mticker.MaxNLocator(3, prune='both'),
        minor_locator=mticker.AutoMinorLocator(),
    )


def rescale_data(interval, data, ax, downsample=0, win_fn=np.mean, **kwargs):
    """Stretch and downsample data points for plotting"""

    fig = ax.figure

    w, h = fig.get_size_inches()
    l = len(interval)

    total_pts = fig.get_dpi() * w / (2**downsample) if downsample is not None else l
    stride = int(l//total_pts)

    # Up- or downsample data
    sample_idx = np.linspace(0, len(data)-1, l).astype(int)

    x = np.arange(interval.start, interval.end)
    y = data[sample_idx]

    if downsample > 0:
        x = x[::stride]
        y = [ win_fn(y[i:i+stride]) for i in np.arange(l)[::stride] ]

    assert len(x) == len(y)
    
    return x, y
