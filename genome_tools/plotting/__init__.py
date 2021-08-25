__all__ =["track"]

import colorsys

import numpy as np
import pandas as pd

from matplotlib.projections import register_projection

import matplotlib.axes as maxes
import matplotlib.axis as maxis
import matplotlib.collections as mcollections
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms
import matplotlib.ticker as mticker

from mpl_toolkits.axes_grid1.inset_locator import inset_axes, InsetPosition
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from genome_tools.helpers import open_file

def scale_lightness(rgb, scale_l):
    # convert rgb to hls
    h, l, s = colorsys.rgb_to_hls(*rgb)
    # manipulate h, l, s values and return as rgb
    return colorsys.hls_to_rgb(h, min(1, l * scale_l), s = s)

def colors_from_cmap(length=50, cmap=None, start=None, stop=None):
    """Return color cycle from a given colormap.
    Parameters
    ----------
    length : int
        The number of colors in the cycle. When `length` is large (> ~10), it
        is difficult to distinguish between successive lines because successive
        colors are very similar.
    cmap : str
        Name of a matplotlib colormap (see matplotlib.pyplot.cm).
    start, stop: 0 <= float <= 1
        Limit colormap to this range (start < stop 1). You should limit the
        range of colormaps with light values (assuming a white background).
        Some colors have default start/stop values (see `CMAP_RANGE`).
    Returns
    -------
    colors : list
        List of RGBA colors.
    See Also
    --------
    cycle_cmap
    """
    if cmap is None:
        cmap = 'viridis'
    if isinstance(cmap, str):
        cmap = getattr(plt.cm, cmap)

    crange = (0, 1) #list(CMAP_RANGE.get(cmap.name, (0, 1)))
    if start is not None:
        crange[0] = start
    if stop is not None:
        crange[1] = stop

    assert 0 <= crange[0] <= 1
    assert 0 <= crange[1] <= 1

    idx = np.linspace(crange[0], crange[1], num=length)
    return cmap(idx)

def cycle_cmap(length=50, cmap=None, start=None, stop=None, ax=None):
    """Set default color cycle of matplotlib based on colormap.
    Note that the default color cycle is **not changed** if `ax` parameter
    is set; only the axes's color cycle will be changed.
    Parameters
    ----------
    length : int
        The number of colors in the cycle. When `length` is large (> ~10), it
        is difficult to distinguish between successive lines because successive
        colors are very similar.
    cmap : str
        Name of a matplotlib colormap (see matplotlib.pyplot.cm).
    start, stop: 0 <= float <= 1
        Limit colormap to this range (start < stop 1). You should limit the
        range of colormaps with light values (assuming a white background).
        Some colors have default start/stop values (see `CMAP_RANGE`).
    ax : matplotlib axes
        If ax is not None, then change the axes's color cycle instead of the
        default color cycle.
    See Also
    --------
    colors_from_cmap, color_mapper
    """
    color_cycle = colors_from_cmap(length, cmap, start, stop)
    ax.set_prop_cycle(color=color_cycle)


# ------------------------




# ------------------------

class AxisAnnotation(mtext.Annotation):
    def __init__(self, text, xy, xytext, xycoords, **kwargs):
        super().__init__(text, xy, xytext=xytext, xycoords=xycoords, **kwargs)

    def get_text_window_extent(self, renderer):
        if self.get_bbox_patch():
            return self.get_bbox_patch().get_window_extent(renderer)
        else:
            return mtext.Text.get_window_extent(self)

    def update_positions(self, renderer):
        """
        """
        x1, y1 = self._get_position_xy(renderer)  # Annotated position.
        # generate transformation,
        self.set_transform(self._get_xy_transform(renderer, self.anncoords))

        if self.arrowprops is None:
            return

        bbox = mtext.Text.get_window_extent(self, renderer)

        d = self.arrowprops.copy()

        ms = d.pop("mutation_scale", self.get_size())
        self.arrow_patch.set_mutation_scale(ms)

        shrink = d.pop('shrink', 0.0)

        if "arrowstyle" not in d:
            self.arrow_patch.set_arrowstyle('-')


        # Pick the corner of the text bbox closest to annotated point.
        xpos = [(bbox.x0, 0), ((bbox.x0 + bbox.x1) / 2, 0.5), (bbox.x1, 1)]
        ypos = [(bbox.y0, 0), ((bbox.y0 + bbox.y1) / 2, 0.5), (bbox.y1, 1)]

        xypos = []

        for i in xpos:
            for j in ypos:
                xypos.append((i, j))

        xmin, ymin = min(xypos, key=lambda v: abs(np.hypot(v[1][0] - y1, v[0][0] - x1)))

        relposx, relposy = xmin[1], ymin[1]
        x, y = xmin[0], ymin[0]

        # x, relposx = min(xpos, key=lambda v: abs(v[0] - x1))
        # y, relposy = min(ypos, key=lambda v: abs(v[0] - y1))

        self._arrow_relpos = (relposx, relposy)
        r = np.hypot(y - y1, x - x1)
        shrink_pts = shrink * r / renderer.points_to_pixels(1)
        self.arrow_patch.shrinkA = self.arrow_patch.shrinkB = shrink_pts

        # adjust the starting point of the arrow relative to the textbox.
        # TODO : Rotation needs to be accounted.
        relposx, relposy = self._arrow_relpos
        x0 = bbox.x0 + bbox.width * relposx
        y0 = bbox.y0 + bbox.height * relposy

        # The arrow will be drawn from (x0, y0) to (x1, y1). It will be first
        # clipped by patchA and patchB.  Then it will be shrunk by shrinkA and
        # shrinkB (in points).  If patch A is not set, self.bbox_patch is used.
        self.arrow_patch.set_positions((x0, y0), (x1, y1))

        if "patchA" in d:
            self.arrow_patch.set_patchA(d.pop("patchA"))
        else:
            if self._bbox_patch:
                self.arrow_patch.set_patchA(self._bbox_patch)
            else:
                if self.get_text() == "":
                    self.arrow_patch.set_patchA(None)
                    return
                pad = renderer.points_to_pixels(4)
                r = mpatches.Rectangle(xy=(bbox.x0 - pad / 2, bbox.y0 - pad / 2),
                              width=bbox.width + pad, height=bbox.height + pad,
                              transform=mtransforms.IdentityTransform(), clip_on=False)
                self.arrow_patch.set_patchA(r)

# ------------------------

class AxisAnnotatedRepositionable(object):
    def __init__(self, pad=0):
        self.pad = pad
        self.annotations = []

    def add_annotation(self, text, pos, offset, loc, **kwargs):
        raise NotImplementedError

    def _reposition_axis_annotations(self, renderer):
        
        order = [x for x,y in sorted(enumerate(self.annotations), key=lambda z: z[1].xy[0])]
        
        mid = len(order)//2
        shift_order = order[:mid][::-1] + order[mid:]

        seen_bboxes = []

        for i, idx in enumerate(shift_order):
            current_artist = self.annotations[idx]
            current_bbox = current_artist.get_text_window_extent(renderer)
        
            if len(seen_bboxes) == 0:
                seen_bboxes.append(current_bbox)
                continue
            
            if self.axis_name == 'x':
                slide = (1,0) if i >= mid else (0, 0)
            elif self.axis_name == 'y':
                slide = (1,1) if i >= mid else (0, 1)
            else:
                raise ValueError(f'Axis {self.axis_name} not supported!')

            slide_dir  = int((slide[0] - 0.5)*2)
            furthest = -slide_dir * float('inf')
            shift = 0
            overlaps = False

            for bbox in seen_bboxes:

                if current_bbox.overlaps(bbox):
                    overlaps = True

                bbox_pt = bbox.get_points()[slide]
                current_bbox_pt = current_bbox.get_points()[1-slide[0], slide[1]]

                if slide_dir * bbox_pt > slide_dir * furthest:
                    furthest = bbox_pt
                    shift = slide_dir * (furthest - current_bbox_pt)
            
            if shift and overlaps:
                xy = np.array(current_artist.get_position())
                xy[slide[1]] += (shift + self.pad) * slide_dir *  (72. / self.figure.get_dpi())
                current_artist.set_position(xy)
                current_artist.update_bbox_position_size(renderer)

            seen_bboxes.append(current_artist.get_text_window_extent(renderer))

# ------------------------

class XAxisAnnotated(maxis.XAxis, AxisAnnotatedRepositionable):
    def __init__(self, *args, **kwargs):
        maxis.XAxis.__init__(self, *args, **kwargs)
        AxisAnnotatedRepositionable.__init__(self, pad=5)
  
    def add_annotation(self, text, pos, offset=15, loc='bottom', **kwargs):
        """
        x: data coordinates

        """
        if loc == 'top':
            y = 1
        elif loc == 'bottom':
            y = 0
            offset = -offset
        else:
            raise ValueError('loc argument must be either top or bottom')

        va = kwargs.pop('va', 'bottom' if loc == 'top' else 'top')
        ha = kwargs.pop('ha', 'center')

        arrowprops = kwargs.pop('arrowprops', dict(arrowstyle='-'))

        t = AxisAnnotation(
            text=text,
            xy=(pos, y), 
            xytext=(0, offset), 
            xycoords=self.axes.get_xaxis_transform('tick1'), 
            textcoords='offset points', 
            ha=ha, va=va, 
            arrowprops=arrowprops,
            **kwargs)
        
        t.set_figure(self.axes.figure)
        self.annotations.append(t)

    def get_annotations(self):
        return self.annotations

    def get_children(self):
        return super().get_children() + [*self.annotations]

    def get_tightbbox(self, renderer, *, for_layout_only=False):
        
        bboxes = [a.get_window_extent(renderer) for a in self.annotations]
        parent_bboxes = super().get_tightbbox(renderer, for_layout_only=for_layout_only)

        if parent_bboxes:
            bboxes.append(parent_bboxes)

        if bboxes:
            return mtransforms.Bbox.union(bboxes)
        else:
            return None

    def draw(self, renderer, *args, **kwargs):
        super().draw(renderer, *args, **kwargs)

        for a in self.annotations:
            a.update_positions(renderer)
            a.update_bbox_position_size(renderer)
                
        self._reposition_axis_annotations(renderer)

        renderer.open_group(__name__, gid=self.get_gid())

        for a in self.annotations:
            a.draw(renderer)

        renderer.close_group(__name__)

# ------------------------

class YAxisAnnotated(maxis.YAxis, AxisAnnotatedRepositionable):
    def __init__(self, *args, **kwargs):
        maxis.YAxis.__init__(self, *args, **kwargs)
        AxisAnnotatedRepositionable.__init__(self, pad=5)
  
    def add_annotation(self, text, pos, offset=15, loc='right', **kwargs):
        """
        x: data coordinates

        """
        if loc == 'left':
            x = 0
            offset = -offset
        elif loc == 'right':
            x = 1
        else:
            raise ValueError('loc argument must be either left or right')

        va = kwargs.pop('va', 'center')
        ha = kwargs.pop('ha', 'right' if loc == 'left' else 'left')

        arrowprops = kwargs.pop('arrowprops', dict(arrowstyle='-'))

        t = AxisAnnotation(
            text=text, 
            xy=(x, pos), 
            xytext=(offset, 0), 
            xycoords=self.axes.get_yaxis_transform('tick1'), 
            textcoords='offset points', 
            horizontalalignment=ha, verticalalignment=va, 
            arrowprops=arrowprops,
            **kwargs)
        
        t.set_figure(self.axes.figure)
        self.annotations.append(t)

    def get_annotations(self):
        return self.annotations

    def _get_annotation_bboxes(self, renderer):
            return [a.get_window_extent(renderer) for a in self.annotations if a.get_visible()]

    def get_ticklabel_extents(self, renderer):
        return super().get_ticklabel_extents(renderer)

    def get_children(self):
        return super().get_children() + [*self.annotations]

    def get_tightbbox(self, renderer, *, for_layout_only=False):
        
        bboxes = [a.get_window_extent(renderer) for a in self.annotations]
        parent_bboxes = super().get_tightbbox(renderer, for_layout_only=for_layout_only)

        if parent_bboxes:
            bboxes.append(parent_bboxes)

        if bboxes:
            return mtransforms.Bbox.union(bboxes)
        else:
            return None
    
    def draw(self, renderer, *args, **kwargs):
        super().draw(renderer, *args, **kwargs)

        for a in self.annotations:
            a.update_positions(renderer)
            a.update_bbox_position_size(renderer)

        self._reposition_axis_annotations(renderer)

        for a in self.annotations:
            a.draw(renderer)

        renderer.close_group(__name__)

# ------------------------

class GenomeAxes(plt.Axes):
    name = 'genome'

    def _init_axis(self):
        self.xaxis = XAxisAnnotated(self)
        self.spines.bottom.register_axis(self.xaxis)
        self.spines.top.register_axis(self.xaxis)
        self.yaxis = YAxisAnnotated(self)
        self.spines.left.register_axis(self.yaxis)
        self.spines.right.register_axis(self.yaxis)
        self._update_transScale()

    def set_interval(self, interval):
        self.set_xlim(interval.start, interval.end)


    def signal(self, data, type='overlay', **kwargs):
        raise NotImplementedError

    def segments(self, intervals, **kwargs):
        raise NotImplementedError

    def seqplot(self):
        raise NotImplemented

    def heatmap(self):
        raise NotImplementedError

    def draw(self, renderer):
        super().draw(renderer)


register_projection(GenomeAxes)
GenomeAxesSubplot = maxes.subplot_class_factory(GenomeAxes)
GenomeAxesSubplot.__module__ = GenomeAxes.__module__

# ------------------------

def clear_spines(ax):
    for loc, spine in ax.spines.items():
        spine.set_visible(False)

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

def format_axes_to_interval(ax, interval, xaxis='bottom', yaxis='left', **kwargs):
    """Format axis appearance"""

    ax.set_xlim(interval.start, interval.end)

    if xaxis:
        ax.xaxis.set_tick_params(direction='out')

        ax.xaxis.set(
            major_locator = mticker.MaxNLocator(3, prune = 'both'),
            major_formatter = mticker.FuncFormatter(get_coord_formatter(interval)),
            minor_locator=mticker.AutoMinorLocator(4))

        if xaxis == 'top':
            ax.xaxis.tick_top()
        elif xaxis == 'bottom':
            ax.xaxis.tick_bottom()
    else:
        ax.xaxis.set_ticks([])
        ax.xaxis.set_ticklabels([])

    if yaxis:
        ax.yaxis.set_tick_params(direction='out')
        
        if yaxis == 'left':
            ax.yaxis.tick_left()
        elif yaxis == 'right':
            ax.yaxis.tick_right()
    else:
        ax.yaxis.set_ticks([])
        ax.yaxis.set_ticklabels([])


def transform_pts(interval, data, ax, downsample=0, win_fn=np.mean, **kwargs):
    """
    Stretch and downsample data points for plotting

    Paramters
    ---------
    interval : genomic_interval
        Region to plot
    data : numpy.ndarray
        Data for cooresponing region
    ax : matplotlib.axes.Axes
        Plotting axes
    downsample : float
        Amount of points to downsample
    win_fun: callable
        Function to apply to down/upsampled windows
    """
    fig = ax.figure
    
    w, h = fig.get_size_inches()
    l = len(interval)

    total_pts = (72. / fig.get_dpi()) * w / (2**downsample) if downsample is not None else l
    stride = int(l//total_pts)

    # Up- or downsample data
    idx = np.linspace(0, len(data)-1, l).astype(int)
    x = np.arange(interval.start, interval.end)
    y = data[idx]

    if downsample:
        nx = x[::stride]
        ny = [win_fn(y[i:i+stride]) for i in np.arange(l)[::stride]]
    else:
        nx, ny = x, y

    assert len(nx) == len(ny)
    return nx, ny

def add_scale_bar(ax):
    scale=100
    bar = AnchoredSizeBar(ax.transData,
            scale, 
            label=f"{scale}nt",
            loc='lower right',
            frameon=False,
            pad=0.4,
            borderpad=1,
            label_top=True,
            fontproperties=dict(size='x-small'))

    bar.patch.set_facecolor('w')
    bar.patch.set_edgecolor('none')
    bar.patch.set_alpha(0.7)
    bar.patch.set_visible(True)
    ax.add_artist(bar)

# ------------------------

ideogram_color_lookup = {
                  'gneg': (1., 1., 1.),
                'gpos25': (.6, .6, .6),
                'gpos50': (.4, .4, .4),
                'gpos75': (.2, .2, .2),
               'gpos100': (0., 0., 0.),
                  'acen': (.8, .4, .4),
                  'gvar': (.8, .8, .8),
                 'stalk': (.9, .9, .9),
               }

def read_ideogram(filepath):
    
    xranges={}
    colors={}
    centromeres={}

    f = open_file(filepath)

    last_chrom=None
    xr=[]
    for line in f:
        chrom, start, stop, label, stain = line.strip().split('\t')
        start = int(start)
        stop = int(stop)
        width = stop - start
    
        if stain=="acen":
            centromeres[chrom] = centromeres.get(chrom, []) + [(start, width)]
            continue
    
        xranges[chrom] = xranges.get(chrom, []) + [(start, width)]
        colors[chrom] = colors.get(chrom, []) + [(ideogram_color_lookup[stain])]

    f.close()

    xranges = xranges
    colors = colors
    centromeres = centromeres

    return xranges, colors, centromeres

def ideogram_plot(data, chrom, pos=None, ax=None, **kwargs):

    if not ax:
        ax = plt.gca()
    fig = ax.figure

    try:
        xranges = data[0][chrom]
        colors = data[1][chrom]
        centromeres = data[2][chrom]
    except:
        print("Error: No chromosome named: {}".format(chrom))
        return

    yranges = (0, 0.5)

    coll = mcollections.BrokenBarHCollection(xranges, yranges, facecolors=colors, edgecolors='black', linewidths=0.5)    
    ax.add_collection(coll)

    if pos:
        ax.axvline(pos, color='red', lw=4)
    w = xranges[-1][0]+xranges[-1][1]

    pad=w*0.05

    ax.set_xlim(0-pad, xranges[-1][0]+xranges[-1][1]+pad)
    ax.xaxis.set_visible(False)

    center = yranges[0] + yranges[1]/2.

    x0, y0 = centromeres[0][0], yranges[0]
    x1, y1 = centromeres[0][0], yranges[1]
    x2, y2 = centromeres[0][0]+centromeres[0][1], center
    cent = mpatches.Polygon(np.array([[x0, y0], [x1, y1], [x2, y2]]), closed=True, 
                    fc=ideogram_color_lookup['acen'], ec='black', linewidth=0.5)
    ax.add_patch(cent)

    x0, y0 = centromeres[1][0], center
    x1, y1 = centromeres[1][0]+centromeres[1][1], yranges[1]
    x2, y2 = centromeres[1][0]+centromeres[1][1], yranges[0]

    cent = mpatches.Polygon(np.array([[x0, y0], [x1, y1], [x2, y2]]), closed=True,
                    fc=ideogram_color_lookup['acen'], ec='black', linewidth=0.5)
    ax.add_patch(cent)

    ax.set_yticks([center])
    ax.set_yticklabels([chrom])
    ax.set_ylim(-0.2, 0.7)

    clear_spines(ax)
# ------------------------

def signal_plot(interval, data, ax=None, **kwargs):
    
    if not ax:
        ax = plt.gca()
    fig = ax.get_figure

    xaxis = kwargs.pop('xaxis', 'bottom')
    yaxis = kwargs.pop('yaxis', 'left')
    
    x, y = transform_pts(interval, data, ax)
    ax.fill_between(x, 0, y, step='mid', **kwargs)

    ax.set_ylim(bottom=0)
    
    format_axes_to_interval(ax, interval, xaxis=xaxis, yaxis=yaxis)

    return ax
# ------------------------

def grouped_heatmap_plot(interval, 
                         data, 
                         group_level=None, 
                         sort_col=None, 
                         hspace=0.01, 
                         ax=None, 
                         cb_loc='bottom', 
                         cb_label='',
                         **kwargs):

    if not ax:
        ax = plt.gca()
    fig = ax.figure

    xaxis = kwargs.pop('xaxis', 'bottom')
    yaxis = kwargs.pop('yaxis', 'left')

    if sort_col is not None:
        df = data.sort_values(by=sort_col)
    else:
        df = data.copy()

    group_sizes = data.groupby(level=group_level).size()

    height_ratios = group_sizes / np.sum(group_sizes)
    total = np.sum(height_ratios) + (hspace * (len(group_sizes)-1)) + 0.15

    height_ratios /= total

    cb_height = 0.15/total
    hspace /= total

    ax.set_xlim(interval.start, interval.end)
    ax.set_ylim(0, 1)

    if cb_loc == 'top':
        axes_pos = 0
        cb_pos = 1.0-cb_height
    else:
        axes_pos = cb_height
        cb_pos = 0

    ticks = []
    ticklabels = []

    for group, group_df in data.groupby(level=group_level):
        
        inset_ax = fig.add_axes([0, 0, 1, 1])
        
        ip = InsetPosition(ax, [0, axes_pos, 1, height_ratios[group]]) #posx, posy, width, height
        inset_ax.set_axes_locator(ip)
        
        hm = inset_ax.pcolormesh(group_df)

        inset_ax.xaxis.set_visible(False)
        inset_ax.yaxis.set_visible(False)    

        mid = axes_pos + (height_ratios[group]/2)
        ticks.append(mid)
        ticklabels.append(group)

        axes_pos += height_ratios[group]+hspace

    # reserve a portion of ax for colorbar
    cb_bbox_ax = fig.add_axes([0, 0, 1, 1], frameon=False)    
    ip = InsetPosition(ax, [0, cb_pos, 0.15, cb_height]) #posx, posy, width, height
    cb_bbox_ax.set_axes_locator(ip)
    cb_bbox_ax.xaxis.set_visible(False)
    cb_bbox_ax.yaxis.set_visible(False)  
    
    cb_ax = inset_axes(cb_bbox_ax, width='100%', height='30%', borderpad=0, loc='center left')
    cb = fig.colorbar(hm, cax=cb_ax, orientation = 'horizontal')
    cb_ax.tick_params(labelsize='xx-small', direction='in')
    cb_ax.yaxis.set_label_position("right")
    cb_ax.set_ylabel(cb_label, fontsize='x-small', rotation=0, ha='left', va='center_baseline')
    

    if hasattr(ax.yaxis, 'add_annotation'):
        [ax.yaxis.add_annotation(text, pos) for pos, text in zip(ticks, ticklabels)]
        ax.yaxis.set_ticks([])
    else:
        ax.yaxis.set_ticks(ticks)
        ax.yaxis.set_ticklabels(ticklabels)

    format_axes_to_interval(ax, interval, xaxis=xaxis, yaxis=yaxis)
    clear_spines(ax)

    return ax

# ------------------------

class row_element:
    def __init__(self, interval, prev=None, next=None, pad=0):
        self.interval = interval
        self.start = interval.start - pad
        self.end = interval.end + pad
        self.prev = prev
        self.next = next

class row:
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

def pack_rows(intervals, pad=5):
    rows = []
    row_idx = []
    curr_row = -1

    for interval in intervals:
        
        e = row_element(interval, pad=pad)

        placed = False 

        for r in rows:
            if r.add(e):
                row_idx.append(r.i)
                placed = True
                break
        
        if placed:
            continue

        curr_row += 1
        r = row(curr_row)
        r.add(e)
        rows.append(r)
        row_idx.append(r.i)

    assert len(row_idx) == len(intervals)
    return row_idx


def segment_plot(interval, segments, pad_points=1, ax=None, **kwargs):
    
    if not ax:
        ax = plt.gca()
    fig = ax.figure

    xaxis = kwargs.pop('xaxis', 'bottom')
    yaxis = kwargs.pop('yaxis', None)

    rectprops = {}
    rectprops['color'] = 'k'
    rectprops['edgecolor'] = 'none'
    # rectprops.update(kwargs)

    ax.set_xlim(interval.start, interval.end)

    trans_left = ax.transData + mtransforms.ScaledTranslation(-1 * pad_points / 72., 0,  ax.figure.dpi_scale_trans)
    trans_right = ax.transData + mtransforms.ScaledTranslation(pad_points / 72., 0,  ax.figure.dpi_scale_trans)

    x0 = ax.transData.inverted().transform(trans_left.transform((interval.start, 0)))[0]
    x1 = ax.transData.inverted().transform(trans_right.transform((interval.start, 0)))[0]

    pad_bp = (x1-x0)//2
    row_idxs = pack_rows(segments, pad=pad_bp)

    patches = []
    for i, row_idx in zip(segments, row_idxs):
        patches.append(mpatches.Rectangle((i.start, row_idx+0.3), i.end-i.start, 0.4, **rectprops))

    pc = mcollections.PatchCollection(patches, match_original=True)
    ax.add_collection(pc)

    ax.set_ylim(min(row_idxs), max(row_idxs)+1)
    ax.invert_yaxis()

    #ax.yaxis.set_ticks(np.array(list(set(row_idxs)))+0.5)
    #ax.yaxis.set_ticklabels([])

    format_axes_to_interval(ax, interval, xaxis=xaxis, yaxis=yaxis)

    return ax

# ------------------------


from .sequence import VOCABS, letter_polygons
from .sequence import add_letter_to_axis

def seq_plot(interval, letter_heights, ax=None, vocab='DNA', **kwargs):

    if not ax:
        ax = plt.gca()
    fig = ax.figure

    xaxis = kwargs.pop('xaxis', 'bottom')
    yaxis = kwargs.pop('yaxis', 'left')

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
        x_pos += interval.start
        for height, letter in letters_and_heights:
            color = VOCABS[vocab][letter]
            polygons = letter_polygons[letter]
            if height > 0:
                add_letter_to_axis(ax, polygons, color, 0.5 + x_pos, y_pos_pos, height)
                y_pos_pos += height
            else:
                add_letter_to_axis(ax, polygons, color, 0.5 + x_pos, y_neg_pos, height)
                y_neg_pos += height
    
    ax.set_ylim(bottom=np.min(letter_heights), top=np.max(letter_heights))
    
    format_axes_to_interval(ax, interval, xaxis=xaxis, yaxis=yaxis)
