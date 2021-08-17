__all__ =["track"]

import numpy as np
import matplotlib.pyplot as plt
import colorsys

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