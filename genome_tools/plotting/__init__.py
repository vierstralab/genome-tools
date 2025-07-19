import numpy as np
import pandas as pd

from collections.abc import Iterable, Sequence, Mapping

from matplotlib.projections import register_projection

import matplotlib.axes as maxes
import matplotlib.axis as maxis
import matplotlib.collections as mcollections
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import matplotlib.transforms as mtransforms
import matplotlib.ticker as mticker

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

# registers custom colormaps
# from .colors import cm
from .utils import rescale_data, format_axes_to_interval, clear_spines, pack_rows

from .gene_annotation import gene_annotation_plot


class AxisAnnotation(mtext.Annotation):
    def __init__(self, text, xy, xytext, xycoords, **kwargs):
        super().__init__(text, xy, xytext=xytext, xycoords=xycoords, **kwargs)

    def get_text_window_extent(self, renderer):
        if self.get_bbox_patch():
            return self.get_bbox_patch().get_window_extent(renderer)
        else:
            return mtext.Text.get_window_extent(self)

    def update_positions(self, renderer):
        """ """
        x1, y1 = self._get_position_xy(renderer)  # Annotated position.
        # generate transformation,
        self.set_transform(self._get_xy_transform(renderer, self.anncoords))

        if self.arrowprops is None:
            return

        bbox = mtext.Text.get_window_extent(self, renderer)

        d = self.arrowprops.copy()

        ms = d.pop("mutation_scale", self.get_size())
        self.arrow_patch.set_mutation_scale(ms)

        shrink = d.pop("shrink", 0.0)

        if "arrowstyle" not in d:
            self.arrow_patch.set_arrowstyle("-")

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
                r = mpatches.Rectangle(
                    xy=(bbox.x0 - pad / 2, bbox.y0 - pad / 2),
                    width=bbox.width + pad,
                    height=bbox.height + pad,
                    transform=mtransforms.IdentityTransform(),
                    clip_on=False,
                )
                self.arrow_patch.set_patchA(r)


# ------------------------


class AxisAnnotatedRepositionable(object):
    def __init__(self, pad=0):
        self.pad = pad
        self.annotations = []

    def add_annotation(self, text, pos, offset, loc, **kwargs):
        raise NotImplementedError

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

    def _get_annotation_bboxes(self, renderer):
        return [
            a.get_window_extent(renderer)
            for a in self.annotations
            if a.get_visible()
        ]

    def _reposition_axis_annotations(self, renderer):
        order = [
            x for x, y in sorted(enumerate(self.annotations), key=lambda z: z[1].xy[0] if self.axis_name == 'x' else z[1].xy[1])
        ]

        mid = len(order) // 2
        shift_order = order[:mid][::-1] + order[mid:]

        seen_bboxes = []

        for i, idx in enumerate(shift_order):
            current_artist = self.annotations[idx]
            current_bbox = current_artist.get_text_window_extent(renderer)

            if len(seen_bboxes) == 0:
                seen_bboxes.append(current_bbox)
                continue

            if self.axis_name == "x":
                slide = (1, 0) if i >= mid else (0, 0)
            elif self.axis_name == "y":
                slide = (1, 1) if i >= mid else (0, 1)
            else:
                raise ValueError(f"Axis {self.axis_name} not supported!")

            slide_dir = int((slide[0] - 0.5) * 2)
            furthest = -slide_dir * float("inf")
            shift = 0
            overlaps = False

            for bbox in seen_bboxes:
                if current_bbox.overlaps(bbox):
                    overlaps = True

                bbox_pt = bbox.get_points()[slide]
                current_bbox_pt = current_bbox.get_points()[1 - slide[0], slide[1]]

                if slide_dir * bbox_pt > slide_dir * furthest:
                    furthest = bbox_pt
                    shift = slide_dir * (furthest - current_bbox_pt)

            if shift and overlaps:
                xy = np.array(current_artist.get_position())
                xy[slide[1]] += (
                    (shift + self.pad) * slide_dir * (72.0 / self.figure.get_dpi())
                )
                current_artist.set_position(xy)
                current_artist.update_bbox_position_size(renderer)

            seen_bboxes.append(current_artist.get_text_window_extent(renderer))

    def draw(self, renderer, *args, **kwargs):
        #super().draw(renderer, *args, **kwargs)

        for a in self.annotations:
            a.update_positions(renderer)
            a.update_bbox_position_size(renderer)

        self._reposition_axis_annotations(renderer)

        renderer.open_group(__name__, gid=self.get_gid())

        for a in self.annotations:
            a.draw(renderer)

        renderer.close_group(__name__)


# ------------------------


class XAxisAnnotated(maxis.XAxis, AxisAnnotatedRepositionable):
    def __init__(self, *args, **kwargs):
        maxis.XAxis.__init__(self, *args, **kwargs)
        AxisAnnotatedRepositionable.__init__(self, pad=5)

    def add_annotation(self, text, pos, offset=15, loc="bottom", **kwargs):
        """
        x: data coordinates

        """
        if loc == "top":
            y = 1
        elif loc == "bottom":
            y = 0
            offset = -offset
        else:
            raise ValueError("`loc` argument must be either top or bottom")

        va = kwargs.pop("va", "bottom" if loc == "top" else "top")
        ha = kwargs.pop("ha", "center")

        arrowprops = kwargs.pop("arrowprops", dict(arrowstyle="-"))

        t = AxisAnnotation(
            text=text,
            xy=(pos, y),
            xytext=(0, offset),
            xycoords=self.axes.get_xaxis_transform("tick1"),
            textcoords="offset points",
            ha=ha,
            va=va,
            arrowprops=arrowprops,
            **kwargs,
        )

        t.set_figure(self.axes.figure)
        self.annotations.append(t)

    def draw(self, renderer, *args, **kwargs):
        AxisAnnotatedRepositionable.draw(self, renderer, *args, **kwargs)

# ------------------------


class YAxisAnnotated(maxis.YAxis, AxisAnnotatedRepositionable):
    def __init__(self, *args, **kwargs):
        maxis.YAxis.__init__(self, *args, **kwargs)
        AxisAnnotatedRepositionable.__init__(self, pad=5)

    def add_annotation(self, text, pos, offset=15, loc="right", **kwargs):
        """
        x: data coordinates

        """
        if loc == "left":
            x = 0
            offset = -offset
        elif loc == "right":
            x = 1
        else:
            raise ValueError("`loc` argument must be either left or right")

        va = kwargs.pop("va", "center")
        ha = kwargs.pop("ha", "right" if loc == "left" else "left")

        arrowprops = kwargs.pop("arrowprops", dict(arrowstyle="-"))

        t = AxisAnnotation(
            text=text,
            xy=(x, pos),
            xytext=(offset, 0),
            xycoords=self.axes.get_yaxis_transform("tick1"),
            textcoords="offset points",
            horizontalalignment=ha,
            verticalalignment=va,
            arrowprops=arrowprops,
            **kwargs,
        )

        t.set_figure(self.axes.figure)
        self.annotations.append(t)

    def draw(self, renderer, *args, **kwargs):
        AxisAnnotatedRepositionable.draw(self, renderer, *args, **kwargs)
# ------------------------


class GenomeAxes(plt.Axes):
    name = "genome"

    def _init_axis(self):
        self.xaxis = XAxisAnnotated(self)
        self.spines.bottom.register_axis(self.xaxis)
        self.spines.top.register_axis(self.xaxis)
        self.yaxis = YAxisAnnotated(self)
        self.spines.left.register_axis(self.yaxis)
        self.spines.right.register_axis(self.yaxis)
        self._update_transScale()


register_projection(GenomeAxes)
GenomeAxesSubplot = maxes.subplot_class_factory(GenomeAxes)
GenomeAxesSubplot.__module__ = GenomeAxes.__module__

# ------------------------

def add_scale_bar(ax):
    scale = 100
    bar = AnchoredSizeBar(
        ax.transData,
        scale,
        label=f"{scale}nt",
        loc="lower right",
        frameon=False,
        pad=0.4,
        borderpad=1,
        label_top=True,
        fontproperties=dict(size="x-small"),
    )

    bar.patch.set_facecolor("w")
    bar.patch.set_edgecolor("none")
    bar.patch.set_alpha(0.7)
    bar.patch.set_visible(True)
    ax.add_artist(bar)


# ------------------------


def signal_plot(interval, data, ax=None, **kwargs):
    if not ax:
        ax = plt.gca()
    fig = ax.get_figure()

    x, y = rescale_data(interval, data, ax, downsample=kwargs.pop("downsample", 0))

    ax.fill_between(x, 0, y, step="mid", **kwargs)

    ax.set_ylim(bottom=0)

    format_axes_to_interval(ax, interval)

    return ax


# ------------------------


def grouped_heatmap_plot(
    interval,
    data,
    group_level=None,
    sort_col=None,
    hspace=0.01,
    cmap="wolfgang_extra",
    ax=None,
    cb_loc="bottom",
    cb_label="",
    **kwargs,
):
    try:
        from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
    except ImportError:
        raise ImportError(
            "Please install mpl version <= 3.10 to use grouped_heatmap_plot."
        )
    if not ax:
        ax = plt.gca()
    fig = ax.get_figure()

    if sort_col is not None:
        df = data.sort_values(by=sort_col)
    else:
        df = data.copy()

    group_sizes = data.groupby(level=group_level).size()

    height_ratios = group_sizes / np.sum(group_sizes)
    total = np.sum(height_ratios) + (hspace * (len(group_sizes) - 1)) + 0.15

    height_ratios /= total

    cb_height = 0.15 / total
    hspace /= total

    ax.set_xlim(interval.start, interval.end)
    ax.set_ylim(0, 1)

    if cb_loc == "top":
        axes_pos = 0
        cb_pos = 1.0 - cb_height
    else:
        axes_pos = cb_height
        cb_pos = 0

    ticks = []
    ticklabels = []

    for group, group_df in data.groupby(level=group_level):
        inset_ax = fig.add_axes([0, 0, 1, 1])

        ip = InsetPosition(
            ax, [0, axes_pos, 1, height_ratios[group]]
        )  # posx, posy, width, height
        inset_ax.set_axes_locator(ip)

        hm = inset_ax.pcolormesh(group_df, cmap=cmap)

        inset_ax.xaxis.set_visible(False)
        inset_ax.yaxis.set_visible(False)

        mid = axes_pos + (height_ratios[group] / 2)
        ticks.append(mid)
        ticklabels.append(group)

        axes_pos += height_ratios[group] + hspace

    # reserve a portion of ax for colorbar
    cb_bbox_ax = fig.add_axes([0, 0, 1, 1], frameon=False)
    ip = InsetPosition(ax, [0, cb_pos, 0.15, cb_height])  # posx, posy, width, height
    cb_bbox_ax.set_axes_locator(ip)
    cb_bbox_ax.xaxis.set_visible(False)
    cb_bbox_ax.yaxis.set_visible(False)

    cb_ax = inset_axes(
        cb_bbox_ax, width="100%", height="30%", borderpad=0, loc="center left"
    )
    cb = fig.colorbar(hm, cax=cb_ax, orientation="horizontal")
    cb_ax.tick_params(labelsize="xx-small", direction="in")
    cb_ax.yaxis.set_label_position("right")
    cb_ax.set_ylabel(
        cb_label, fontsize="x-small", rotation=0, ha="left", va="center_baseline"
    )

    if hasattr(ax.yaxis, "add_annotation"):
        [ax.yaxis.add_annotation(text, pos) for pos, text in zip(ticks, ticklabels)]
        ax.yaxis.set_ticks([])
    else:
        ax.yaxis.set_ticks(ticks)
        ax.yaxis.set_ticklabels(ticklabels)

    format_axes_to_interval(ax, interval, xaxis=xaxis, yaxis=yaxis)
    clear_spines(ax)

    return ax


# ------------------------


def segment_plot(interval, segments, pad_points=1, ax=None, rect_height=0.4, **kwargs):
    if not ax:
        ax = plt.gca()

    assert 0 < rect_height <= 1
    rect_pad = (1 - rect_height) / 2

    rectprops = {}
    rectprops["color"] = "k"
    rectprops["edgecolor"] = "none"
    rectprops.update(kwargs)

    ax.set_xlim(interval.start, interval.end)

    trans_left = ax.transData + mtransforms.ScaledTranslation(
        -1 * pad_points / 72.0, 0, ax.figure.dpi_scale_trans
    )
    trans_right = ax.transData + mtransforms.ScaledTranslation(
        pad_points / 72.0, 0, ax.figure.dpi_scale_trans
    )

    x0 = ax.transData.inverted().transform(trans_left.transform((interval.start, 0)))[0]
    x1 = ax.transData.inverted().transform(trans_right.transform((interval.start, 0)))[0]

    pad_bp = (x1 - x0) // 2
    row_indices = pack_rows(segments, pad=pad_bp)

    summit_lines = []
    patches = []
    for interval_i, row_index in row_indices.items():
        interval_rectprops = rectprops.copy()
        interval_rectprops.update(getattr(interval_i, "rectprops", {}))
        patches.append(
            mpatches.Rectangle(
                (interval_i.start, row_index + rect_pad), interval_i.end - interval_i.start, rect_height, **interval_rectprops
            )
        )
        if hasattr(interval_i, "summit"):
            summit_lines.append(
                mlines.Line2D(
                    [interval_i.summit, interval_i.summit],
                    [row_index + rect_pad - rect_height*0.2, row_index + rect_pad + rect_height*1.2],
                    color="k",
                    lw=0.25,
                )
        )

    pc = mcollections.PatchCollection(patches, match_original=True)
    ax.add_collection(pc)

    for line in summit_lines:
        ax.add_line(line)

    ax.set_ylim(
        min(row_indices.values(), default=0),
        max(row_indices.values(), default=1) + 1
    )
    ax.invert_yaxis()

    # ax.yaxis.set_ticks(np.array(list(set(row_idxs)))+0.5)
    # ax.yaxis.set_ticklabels([])

    format_axes_to_interval(ax, interval)

    return ax


# ------------------------
