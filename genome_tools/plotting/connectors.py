
from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector, BboxConnectorPatch

# Connect two bounding boxes in possibly different axes
def connect_bbox(bbox1, bbox2, loc1a, loc2a, loc1b, loc2b, prop_lines, prop_patches=None):
    if prop_patches is None:
        prop_patches = prop_lines.copy()
        prop_patches["alpha"] = prop_patches.get("alpha", 1)*0.2

    c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
    c1.set_clip_on(False)

    c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
    c2.set_clip_on(False)

    bbox1_patch= BboxPatch(bbox1, **prop_patches)
    bbox2_patch = BboxPatch(bbox2, **prop_patches)

    # connects two bbox-es with a quadrilateral
    p = BboxConnectorPatch(bbox1, bbox2, loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b, **prop_patches)
    p.set_fill(True)
    p.set_clip_on(False)

    return c1, c2, bbox1_patch, bbox2_patch, p

def zoom_effect(ax_zoomed, ax_origin, xlims = None, orientation='below', **kwargs):
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

        bbox_zoomed=ax_zoomed.bbox
        bbox_origin=TransformedBbox(ax_zoomed.viewLim, transform)
    else:
        transform_zoomed=blended_transform_factory(ax_zoomed.transData, ax_zoomed.transAxes)
        transform_origin=blended_transform_factory(ax_origin.transData, ax_origin.transAxes)
    
        bbox_zoomed=TransformedBbox(Bbox.from_extents(xlims[0], 0, xlims[1], 1), transform_zoomed)
        bbox_origin=TransformedBbox(Bbox.from_extents(xlims[0], 0, xlims[1], 1), transform_origin)

    prop_patches = kwargs.copy()
    prop_patches["ec"] = "none"
    prop_patches["alpha"] = 0.2

    if orientation=='below':
        loc1a=2
        loc2a=3
        loc1b=1
        loc2b=4
    elif orientation=='above':
        loc1a=3
        loc2a=2
        loc1b=4
        loc2b=1
    else:
        raise Exception("orientation '%s' not recognized" % orientation)

    c1, c2, bbox_zoomed_patch, bbox_origin_patch, p = \
        connect_bbox(bbox_zoomed, bbox_origin,
                     loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                     prop_lines=kwargs, prop_patches=prop_patches)

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


    c1, c2, bbox_patch1, bbox_patch2, p = \
        connect_bbox(mybbox1, mybbox2,
                     loc1a=2, loc2a=3, loc1b=1, loc2b=4,
                     prop_lines=kwargs, prop_patches=prop_patches)

    #ax1.add_patch(bbox_patch1)
    #ax2.add_patch(bbox_patch2)
    ax2.add_patch(c1)
    ax2.add_patch(c2)
    ax2.add_patch(p)

    return c1, c2, bbox_patch1, bbox_patch2, p

