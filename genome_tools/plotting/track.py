# Copyright 2016 Jeff Vierstra

import matplotlib.pyplot as plt
import matplotlib.axes as maxes
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches

from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from genome_tools.plotting import scale_lightness, cycle_cmap
from genome_tools.plotting.annotation import annotate, annotate_span
from genome_tools.plotting.sequence import add_letter_to_axis, letter_polygons, VOCABS

import numpy as np


def get_coord_formatter(interval):
    if len(interval) > 1e5:
        fmt='{:0.2f} mb'
        denom=1e6
    elif len(interval) > 1e3:
        fmt='{:0.2f} kb'
        denom=1e3
    else:
        fmt='{:,.0f}'
        denom=1

    def tick_format_func(x, pos):
        return fmt.format(x/denom) 
    
    return tick_format_func

class base_track(object):
    def __init__(self, interval, ax=None):
        if ax is None:
            self.ax = plt.gca()
        else:
            self.ax = ax

        self.interval = interval
        self.data = []
        self.data_kwargs = []

    def transform_pts(self, data, downsample=0, win_fn=np.mean, **kwargs):
        """"""
        fig = self.ax.get_figure()
        
        w, h = fig.get_size_inches()
        dpi = fig.get_dpi()
        N = len(self.interval)

        total_pts = dpi*w / (2**downsample) if downsample is not None else N
        stride = int(N//total_pts)
 
        # Up- or downsample data
        idx = np.linspace(0, len(data)-1, N).astype(int)
        x = np.arange(self.interval.start, self.interval.end)
        y = data[idx]

        if downsample:
            nx = x[::stride]
            ny = [win_fn(y[i:i+stride]) for i in np.arange(N)[::stride]]
        else:
            nx, ny = x, y

        assert len(nx) == len(ny)

        return nx, ny

    def format_spines(self, remove=['top', 'right']):
        """Remove spines"""
        for spine in remove:
            self.ax.spines[spine].set_color('none')
            self.ax.spines[spine].set_visible('none')

    def format_axes(self):
        """Format axis appearance"""
        self.ax.xaxis.set_tick_params(direction='out')
        self.ax.yaxis.set_tick_params(direction='out')

        self.ax.xaxis.set(
            major_locator = mticker.MaxNLocator(3, prune = 'both'),
            major_formatter = mticker.FuncFormatter(get_coord_formatter(self.interval)),
            minor_locator=mticker.AutoMinorLocator(4))

        self.ax.set_xlim(self.interval.start, self.interval.end)
    
    def show_grid(self, axis='both'):
        self.ax.grid(axis)
        return self

    def show_legend(self, **kwargs):
        self.ax.legend(**kwargs)
        return self
    
    def show_scale_bar(self, **kwargs):
        scale=100
        bar = AnchoredSizeBar(self.ax.transData,
                scale, 
                label=f"{scale}nt",
                loc=2,
                frameon=False)
        self.ax.add_artist(bar)
        return self

    def add_annotation(self, data_xy, labels, **kwargs):
        annotate(self.ax, labels, data_xy, **kwargs)
        return self

    def add_annotation_span(self, intervals, **kwargs):
        annotate_span(self.ax, intervals, **kwargs)
        return self

    def render(self, **kwargs):
        self.format_spines()
        self.format_axes()

        self.ax.patch.set_color('none')

        if 'ylim' in kwargs:
            self.ax.set_ylim(kwargs.pop('ylim'))

        return self

class genome_track(base_track):
    """
    """
    def __init__(self, interval, ax=None):
        super(genome_track, self).__init__(interval, ax=ax)

    def add_data(self, dataset, **kwargs):
        self.data.append(dataset)
        self.data_kwargs.append(kwargs)
        return self

    def render(self, **kwargs):
        super(genome_track, self).render(**kwargs)

        for d, k in zip(self.data[::-1], self.data_kwargs[::-1]):
            x, y = self.transform_pts(d, **kwargs)
        
            
            fb = self.ax.fill_between(x, 0, y, step='mid', **k)
            
            #darken_col = scale_lightness(fb.get_facecolor()[0][:-1], 0.4)
            #self.ax.step(x, y, color=darken_col, lw=0.5)

        return self

class segment_track(base_track):
    def __init__(self, interval, ax=None, **kwargs):
        super(segment_track, self).__init__(interval, ax=ax, **kwargs)

    def add_data(self, dataset, **kwargs):
        self.data.append(dataset)
        self.data_kwargs.append(kwargs)
        return self

    def render(self, **kwargs):
        super(segment_track, self).render()

        if 'cmap' in kwargs:
            cmap = kwargs.pop('cmap', None)
            cycle_cmap(len(self.data), cmap=cmap, ax=self.ax)

        prop_cycler = self.ax._get_lines.prop_cycler

        for i, segments in enumerate(self.data[::-1]):
            self.ax.broken_barh(
                [(s.start, len(s)) for s in segments['data']],
                (i-0.4, 0.8), fc=next(prop_cycler)['color'],
                **kwargs)
        
        self.ax.yaxis.set_ticks(range(len(self.data)))
        self.ax.yaxis.set_ticklabels([s['name'] for s in self.data[::-1]])
        for label in self.ax.get_yticklabels():
            label.set_fontweight('bold')
                
        return self

class heatmap_track(base_track):
    def __init__(self, interval, data, ax=None, **kwargs):
        super(heatmap_track, self).__init__(interval, data, ax=ax, **kwargs)

    def render(self):
        pass

class seq_track(base_track):
    def __init__(self, interval, ax=None):
        super(seq_track, self).__init__(interval, ax=ax)

    def render(self, letter_heights, vocab='DNA'):
        
        assert letter_heights.shape[1] == len(VOCABS[vocab])
        x_range = [1, letter_heights.shape[0]]
        pos_heights = np.copy(letter_heights)
        pos_heights[letter_heights < 0] = 0
        neg_heights = np.copy(letter_heights)
        neg_heights[letter_heights > 0] = 0

        for x_pos, heights in enumerate(letter_heights):
            letters_and_heights = sorted(zip(heights, list(VOCABS[vocab].keys())))
            y_pos_pos = 0.0
            y_neg_pos = 0.0
            for height, letter in letters_and_heights:
                color = VOCABS[vocab][letter]
                polygons = letter_polygons[letter]
                if height > 0:
                    add_letter_to_axis(self.ax, polygons, color, 0.5 + x_pos, y_pos_pos, height)
                    y_pos_pos += height
                else:
                    add_letter_to_axis(self.ax, polygons, color, 0.5 + x_pos, y_neg_pos, height)
                    y_neg_pos += height

        # if add_hline:
        #     ax.axhline(color="black", linewidth=1)
        self.ax.set_xlim(x_range[0] - 1, x_range[1] + 1)
        self.ax.grid(False)
        self.ax.set_xticks(list(range(*x_range)) + [x_range[-1]])
        self.ax.set_aspect(aspect='auto', adjustable='box')
        self.ax.autoscale_view()

        