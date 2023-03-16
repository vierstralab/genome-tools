import numpy as np

import matplotlib.pyplot as plt
import matplotlib.collections as mcollections
import matplotlib.patches as mpatches

from .utils import clear_spines

from genome_tools.plotting.colors.cm import get_vocab_color
from genome_tools.helpers import open_file


def read_ideogram(filepath):
    xranges = {}
    colors = {}
    centromeres = {}

    with open_file(filepath) as f:
        last_chrom = None
        xr = []
        for line in f:
            chrom, start, stop, label, stain = line.strip().split("\t")
            start = int(start)
            stop = int(stop)
            width = stop - start

            if stain == "acen":
                centromeres[chrom] = centromeres.get(chrom, []) + [(start, width)]
                continue

            xranges[chrom] = xranges.get(chrom, []) + [(start, width)]
            colors[chrom] = colors.get(chrom, []) + [
                (get_vocab_color(stain, "ideogram"))
            ]

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

    coll = mcollections.BrokenBarHCollection(
        xranges, yranges, facecolors=colors, edgecolors="black", linewidths=0.5
    )
    ax.add_collection(coll)

    if pos:
        ax.axvline(pos, color="red", lw=4)
    w = xranges[-1][0] + xranges[-1][1]

    pad = w * 0.05

    ax.set_xlim(0 - pad, xranges[-1][0] + xranges[-1][1] + pad)
    ax.xaxis.set_visible(False)

    center = yranges[0] + yranges[1] / 2.0

    x0, y0 = centromeres[0][0], yranges[0]
    x1, y1 = centromeres[0][0], yranges[1]
    x2, y2 = centromeres[0][0] + centromeres[0][1], center
    cent = mpatches.Polygon(
        np.array([[x0, y0], [x1, y1], [x2, y2]]),
        closed=True,
        fc=get_vocab_color("acen", "ideogram"),
        ec="black",
        linewidth=0.5,
    )
    ax.add_patch(cent)

    x0, y0 = centromeres[1][0], center
    x1, y1 = centromeres[1][0] + centromeres[1][1], yranges[1]
    x2, y2 = centromeres[1][0] + centromeres[1][1], yranges[0]

    cent = mpatches.Polygon(
        np.array([[x0, y0], [x1, y1], [x2, y2]]),
        closed=True,
        fc=get_vocab_color("acen", "ideogram"),
        ec="black",
        linewidth=0.5,
    )
    ax.add_patch(cent)

    ax.set_yticks([center])
    ax.set_yticklabels([chrom])
    ax.set_ylim(-0.2, 0.7)

    clear_spines(ax)
