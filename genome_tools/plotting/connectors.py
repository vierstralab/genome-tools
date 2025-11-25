from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import (
    BboxPatch,
    BboxConnector,
    BboxConnectorPatch,
)

from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import numpy as np


# Connect two bounding boxes in possibly different axes
def connect_bbox(
    bbox1, bbox2, loc1a, loc2a, loc1b, loc2b, prop_lines, prop_patches=None
):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1) * 0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)

    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox1_patch = BboxPatch(bbox1, **prop_patches)
    bbox2_patch = BboxPatch(bbox2, **prop_patches)

    # connects two bbox-es with a quadrilateral
    p = BboxConnectorPatch(
        bbox1, bbox2, loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b, **prop_patches
    )
    p.set_fill(True)
    p.set_clip_on(False)

    return c1, c2, bbox1_patch, bbox2_patch, p


def zoom_effect(ax_zoomed, ax_origin, xlims=None, orientation="below", **kwargs):
    """
    ax_zoomed : zoomed axes
    ax_origin:  the main axes
    (xmin,xmax) : the limits of the colored area in both plot axes.

    connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
    be marked.  The keywords parameters will be used ti create
    patches.

    """
    if xlims is None:
        tt = ax_zoomed.transScale + (ax_zoomed.transLimits + ax_origin.transAxes)
        transform = blended_transform_factory(ax_origin.transData, tt)

        bbox_zoomed = ax_zoomed.bbox
        bbox_origin = TransformedBbox(ax_zoomed.viewLim, transform)
    else:
        transform_zoomed = blended_transform_factory(
            ax_zoomed.transData, ax_zoomed.transAxes
        )
        transform_origin = blended_transform_factory(
            ax_origin.transData, ax_origin.transAxes
        )

        bbox_zoomed = TransformedBbox(
            Bbox.from_extents(xlims[0], 0, xlims[1], 1), transform_zoomed
        )
        bbox_origin = TransformedBbox(
            Bbox.from_extents(xlims[0], 0, xlims[1], 1), transform_origin
        )

    prop_patches = kwargs.copy()
    prop_patches["ec"] = "none"
    prop_patches["alpha"] = 0.2

    if orientation == "below":
        loc1a = 2
        loc2a = 3
        loc1b = 1
        loc2b = 4
    elif orientation == "above":
        loc1a = 3
        loc2a = 2
        loc1b = 4
        loc2b = 1
    else:
        raise Exception("orientation '%s' not recognized" % orientation)

    c1, c2, bbox_zoomed_patch, bbox_origin_patch, p = connect_bbox(
        bbox_zoomed,
        bbox_origin,
        loc1a=loc1a,
        loc2a=loc2a,
        loc1b=loc1b,
        loc2b=loc2b,
        prop_lines=kwargs,
        prop_patches=prop_patches,
    )

    ax_zoomed.add_patch(bbox_zoomed_patch)
    ax_origin.add_patch(bbox_origin_patch)
    ax_origin.add_patch(c1)
    ax_origin.add_patch(c2)
    ax_origin.add_patch(p)

    return c1, c2, bbox_zoomed_patch, bbox_origin_patch, p


def zoom_effect02(ax1, ax2, **kwargs):
    """
    ax1 : the main axes
    ax1 : the zoomed axes

    Similar to zoom_effect01.  The xmin & xmax will be taken from the
    ax1.viewLim.
    """

    tt = ax1.transScale + (ax1.transLimits + ax2.transAxes)
    trans = blended_transform_factory(ax2.transData, tt)

    mybbox1 = ax1.bbox
    mybbox2 = TransformedBbox(ax1.viewLim, trans)

    prop_patches = kwargs.copy()
    prop_patches["ec"] = "none"
    prop_patches["alpha"] = 0.2

    c1, c2, bbox_patch1, bbox_patch2, p = connect_bbox(
        mybbox1,
        mybbox2,
        loc1a=2,
        loc2a=3,
        loc1b=1,
        loc2b=4,
        prop_lines=kwargs,
        prop_patches=prop_patches,
    )

    # ax1.add_patch(bbox_patch1)
    # ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p


def get_fig_coords(ax, x, y, coord_system='data'):
    """
    Transform coordinates (x, y) from a specified coordinate system to figure coordinates.

    Parameters:
        ax: Matplotlib Axes object.
        x, y: Coordinates in the specified coordinate system of the Axes.
        coord_system: The coordinate system of the input points. Options are:
                      - 'data': Data coordinates of the Axes.
                      - 'axes': Axes coordinates (0 to 1 within the Axes).
                      - 'figure': Figure coordinates (0 to 1 within the Figure).
                      - 'display': Display (screen) coordinates.

        x and y can be scalars or array-like.

    Returns:
        Tuple of (x_fig, y_fig) in figure coordinates.
    """
    x = np.atleast_1d(x)
    y = np.atleast_1d(y)
    fig = ax.figure

    if coord_system == 'data':
        coords_disp = ax.transData.transform(np.column_stack([x, y]))
    elif coord_system == 'axes':
        coords_disp = ax.transAxes.transform(np.column_stack([x, y]))
    elif coord_system == 'figure':
        coords_disp = fig.transFigure.transform(np.column_stack([x, y]))
    elif coord_system == 'display':
        coords_disp = np.column_stack([x, y])  # Already in display coordinates
    else:
        raise ValueError(f"Invalid coordinate system: {coord_system}")

    # Transform from display to figure coordinates
    coords_fig = fig.transFigure.inverted().transform(coords_disp)
    x_fig, y_fig = coords_fig[:, 0], coords_fig[:, 1]
    return x_fig, y_fig


def get_axes_vertical_extent_in_fig_coords(ax):
    """
    Get the y positions of the bottom and top of an axis in figure coordinates.

    Parameters:
        ax: Matplotlib Axes object.

    Returns:
        y_bottom_fig, y_top_fig: y-positions in figure coordinates.
    """
    bbox = ax.get_position()
    y_bottom_fig = bbox.y0
    y_top_fig = bbox.y1
    return y_bottom_fig, y_top_fig


def get_connector_points(ax_top, ax_bottom, *X, coord_system='data',
                         extend_to_top=False, extend_to_bottom=False):
    """
    Get the points in figure coordinates
    for connecting two axes vertically
    for each x-coordinate in X.
    Each x-coordinate is connected by a line
    from the top of ax_bottom to the bottom of ax_top.
    If extend_to_top is True, the line extends vertically
    to the top of ax_top. If extend_to_bottom is True,
    the line extends vertically to the bottom of ax_bottom.
    """
    X = np.sort(np.ravel(X))
    
    if ax_top.figure != ax_bottom.figure:
        raise ValueError("ax1 and ax2 must be in the same figure.")

    x_ax1_fig, _ = get_fig_coords(ax_top, X, np.zeros_like(X), coord_system=coord_system)
    x_ax2_fig, _ = get_fig_coords(ax_bottom, X, np.zeros_like(X), coord_system=coord_system)

    y_ax1_bottom_fig, y_ax1_top_fig = get_axes_vertical_extent_in_fig_coords(ax_top)
    y_ax2_bottom_fig, y_ax2_top_fig = get_axes_vertical_extent_in_fig_coords(ax_bottom)

    y_vals = ([y_ax2_bottom_fig] if extend_to_bottom else []) + [y_ax2_top_fig, y_ax1_bottom_fig] + ([y_ax1_top_fig] if extend_to_top else [])
    y_vals = np.tile(y_vals, X.shape[0])

    x_vals = [x_ax2_fig] * (2 if extend_to_bottom else 1) + [x_ax1_fig] * (2 if extend_to_top else 1)
    x_vals = np.ravel(np.column_stack(x_vals))

    shape = 2 + extend_to_bottom + extend_to_top

    return x_vals.reshape((-1, shape)), y_vals.reshape((-1, shape))


def connect_axes_lines(ax_top, ax_bottom, x, coord_system='data', ls='--', extend_to_top=False, extend_to_bottom=False, **kwargs):
    """
    Draw lines connecting corresponding points on two axes vertically.
    """
    points_x, points_y = get_connector_points(ax_top, ax_bottom, *x, coord_system=coord_system,
                                  extend_to_top=extend_to_top, extend_to_bottom=extend_to_bottom)
    
    fig = ax_top.figure

    for px, py in zip(points_x, points_y):
        line = Line2D(px, py, transform=fig.transFigure, ls=ls, **kwargs)
        fig.add_artist(line)


def connect_axes_area(ax_top, ax_bottom, x1, x2, coord_system='data', closed=True, facecolor='grey', alpha=0.3, extend_to_top=False, extend_to_bottom=False, **kwargs):
    """
    Connect two axes with a shaded area between specified x1 and x2 coordinates.
    """
    fig = ax_top.figure

    left,right = np.sort([x1, x2])

    (x_left, x_right), (y_left, y_right) = get_connector_points(ax_top, ax_bottom, left, right, coord_system=coord_system,
                                extend_to_top=extend_to_top, extend_to_bottom=extend_to_bottom)

    x = np.concatenate([x_left, x_right[::-1]])
    y = np.concatenate([y_left, y_right[::-1]])

    vertices = np.column_stack([x, y])

    # Create and add the polygon
    polygon = Polygon(vertices, closed=closed, facecolor=facecolor, alpha=alpha, transform=fig.transFigure, **kwargs)
    fig.patches.append(polygon)