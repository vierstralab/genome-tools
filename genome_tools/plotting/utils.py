import numpy as np
from genome_tools import GenomicInterval
from typing import List
import matplotlib.ticker as mticker

__all__ = [
    "clear_spines",
    "get_genomic_axis_tick_formatter",
    "format_axes_to_interval",
    "rescale_data",
]


def clear_spines(ax):
    """Clear spines from plot"""

    for loc, spine in ax.spines.items():
        spine.set_visible(False)

    ax.patch.set_visible(False)


def get_genomic_axis_tick_formatter(interval):
    """Return a tick formattor for the genomic coordinate axis"""

    if len(interval) > 1e5:
        fmt = "{:0.2f} Mb"
        denom = 1e6
    else:
        fmt = "{:,.0f}"
        denom = 1

    def _func(x, pos):
        return fmt.format(x / denom)

    return mticker.FuncFormatter(_func)


def format_axes_to_interval(ax, interval, axis="both", **kwargs):
    """Format axis appearance"""

    ax.set_xlim(interval.start, interval.end)

    ax.tick_params(axis=axis, direction="out")

    if axis in ("x", "both"):
        ax.xaxis.set(
            major_locator=mticker.MaxNLocator(3, prune="both"),
            minor_locator=mticker.AutoMinorLocator(4),
            major_formatter=get_genomic_axis_tick_formatter(interval),
        )
    if axis in ("y", "both"):
        ax.yaxis.set(
            major_locator=mticker.MaxNLocator(3, prune="both"),
        )
    return ax


def rescale_data(interval, data, ax, downsample=0, win_fn=np.mean, **kwargs):
    """Stretch and downsample data points for plotting"""

    fig = ax.figure

    w, h = fig.get_size_inches()
    l = len(interval)

    total_pts = fig.get_dpi() * w / (2**downsample) if downsample is not None else l
    stride = int(l // total_pts)

    # Up- or downsample data
    sample_idx = np.linspace(0, len(data) - 1, l).astype(int)

    x = np.arange(interval.start, interval.end)
    y = data[sample_idx]

    if downsample > 0:
        x = x[::stride]
        y = [win_fn(y[i : i + stride]) for i in np.arange(l)[::stride]]

    assert len(x) == len(y)

    return x, y


class RowElement:
    def __init__(self, interval, prev=None, next=None, pad=0):
        self.interval = interval
        self.start = interval.start - pad
        self.end = interval.end + pad
        self.prev = prev
        self.next = next


class Row:
    def __init__(self, i):
        self.i = i
        self.first = None
        self.last = None

    def add(self, e):
        if self.first is None:
            e.prev = None
            e.next = None
            self.first = self.last = e
            return 1

        curr = self.first
        while curr and curr.start < e.start:
            curr = curr.next

        if curr is None:
            if self.last.end < e.start:
                e.prev = self.last
                e.next = None
                self.last.next = e
                self.last = e
                return 1
            else:
                return 0

        prev = curr.prev
        if prev is None:
            if e.end < curr.start:
                e.prev = None
                e.next = curr
                curr.prev = e
                self.first = e
                return 1
            else:
                return 0

        if prev.end < e.start and curr.start > e.end:
            e.prev = prev
            e.next = curr
            prev.next = e
            curr.prev = e
            return 1

        return 0


def pack_rows(intervals: List[GenomicInterval], pad=5):
    rows: List[Row] = []
    row_indices = {}
    curr_row = -1

    result_intervals = []
    for interval in intervals:
        e = RowElement(interval, pad=pad)

        placed = False

        for r in rows:
            if r.add(e):
                row_indices[e.interval] = r.i
                placed = True
                break

        if placed:
            continue

        curr_row += 1
        r = Row(curr_row)
        r.add(e)
        rows.append(r)
        row_indices[e.interval] = r.i
        interval_copy = interval.copy()
        interval_copy.row_index = r.i
        result_intervals.append(interval_copy)

    assert len(row_indices) == len(intervals)
    return row_indices, result_intervals
