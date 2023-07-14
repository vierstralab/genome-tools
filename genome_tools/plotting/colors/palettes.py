import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.cm as mcm

from cycler import cycler
from itertools import cycle

import numpy as np

from .cm import COLOR_MAPS

__all__ = ["get_color_cycle", "set_palette", "color_palette"]

QUAL_PALETTES = [
    "zissou",
    "darjeeling",
    "lawhoops",
]
QUAL_PALETTES_SIZES = {name: len(COLOR_MAPS[name]) for name in QUAL_PALETTES}


def get_color_cycle():
    """Returns the current color cycle from matplotlib prop cycler"""
    cyl = mpl.rcParams["axes.prop_cycler"]
    return cyl.by_key()["color"] if "color" in cyl.keys else ["k"]


class _ColorPalette(list):
    """Enables color palette used with `with` command"""

    def __enter__(self):
        self.orig_palette = color_palette()
        set_palette(self)
        return self

    def __exit__(self):
        set_palette(self.orig_palette)


def mpl_cmap_palette(cmap, n_colors, as_cmap=False):
    """Returns N colors (evenly distributed) from a colormap"""

    if isinstance(cmap, str):
        cmap = mcm.get_cmap(cmap)

    bins = np.linspace(0, 1, int(n_colors) + 2)[1:-1]
    colors = list(map(tuple, cmap(bins)[:, :3]))

    if as_cmap:
        return cmap
    else:
        return _ColorPalette(colors)


def discrete_cmap(cmap, n_colors):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    if isinstance(cmap, str):
        cmap = mcm.get_cmap(cmap)

    color_list = mpl_cmap_palette(cmap, n_colors, as_cmap=False)
    return mcolors.ListedColormap(color_list)


def color_palette(palette=None, n_colors=None, as_cmap=False):
    if palette is None:
        palette = get_color_cycle()
        if n_colors is None:
            n_colors = len(palette)

    elif not isinstance(palette, str):
        palette = palette
        if n_colors is None:
            n_colors = len(palette)

    else:
        if n_colors is None:
            n_colors = QUAL_PALETTES_SIZES.get(palette, 6)

        if palette in QUAL_PALETTES:
            palette = COLOR_MAPS[palette]
        else:
            try:
                palette = mpl_cmap_palette(palette, n_colors, as_cmap=as_cmap)
            except:
                raise ValueError(f"{palette} is not a valid palette name.")

    if not as_cmap:
        palette_cycle = cycle(palette)
        palette = [next(palette_cycle) for _ in range(n_colors)]

        try:
            palette = map(mcolors.to_rgb, palette)
            palette = _ColorPalette(palette)
        except ValueError:
            raise ValueError(
                f"Could not generate a palette for {palette}. Check color formats"
            )

    return palette


def set_palette(palette, n_colors):
    """Sets the palette to the default matplotlib prop cycler"""

    colors = color_palette(palette, n_colors)
    cyl = cycler("color", colors)
    mpl.rcParams["axes.prop_cycle"] = cyl
    mpl.rcParams["patch.facecolor"] = colors[0]
