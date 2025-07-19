import matplotlib.colors as mcolors
import matplotlib.cm as mcm
import matplotlib as mpl

from collections.abc import Iterable

__all__ = [
    "COLOR_MAPS",
    "VOCAB_COLOR_MAPS",
    "get_vocab_color",
    "map_vocab_color",
    "discrete_cmap",
]

VOCAB_COLOR_MAPS = dict(
    dna={
        "A": "green",
        "C": "blue",
        "G": "orange",
        "T": "red",
    },
    rna={
        "A": "green",
        "C": "blue",
        "G": "orange",
        "U": "red",
    },
    aa={
        "A": "#CCFF00",
        "B": "orange",
        "C": "#FFFF00",
        "D": "#FF0000",
        "E": "#FF0066",
        "F": "#00FF66",
        "G": "#FF9900",
        "H": "#0066FF",
        "I": "#66FF00",
        "K": "#6600FF",
        "L": "#33FF00",
        "M": "#00FF00",
        "N": "#CC00FF",
        "P": "#FFCC00",
        "Q": "#FF00CC",
        "R": "#0000FF",
        "S": "#FF3300",
        "T": "#FF6600",
        "V": "#99FF00",
        "W": "#00CCFF",
        "Y": "#00FFCC",
        "Z": "blue",
    },
    ideogram={
        "gneg": (1.0, 1.0, 1.0),
        "gpos25": (0.6, 0.6, 0.6),
        "gpos50": (0.4, 0.4, 0.4),
        "gpos75": (0.2, 0.2, 0.2),
        "gpos100": (0.0, 0.0, 0.0),
        "acen": (0.8, 0.4, 0.4),
        "gvar": (0.8, 0.8, 0.8),
        "stalk": (0.9, 0.9, 0.9),
    },
    gene_annotation={
        "exon": "gold",
        "utr": "black",
        "cds": "black",
    },
    # DHS component annotations from Meuleman et al., 2020 Nature
    # meuleman2020 = {
    #     '': '#FFE500',
    #     '': '#FE8102',
    #     '': '#FF0000',
    #     '': '#07AF00',
    #     '': '#4C7D14',
    #     '': '#414613',
    #     '': '#05C1D9',
    #     '': '#0467FD',
    #     '': '#009588',
    #     '': '#BB2DD4',
    #     '': '#7A00FF',
    #     '': '#4A6876',
    #     '': '#08245B',
    #     '': '#B9461D',
    #     '': '#692108',
    #     '': '#C3C3C3',
    # },
)


def map_vocab_color(x, vocab, default="black"):
    """Map vocabulary colors for a list of values"""

    _func = lambda x: get_vocab_color(x, vocab, default)
    if isinstance(x, Iterable):
        return list(map(_func, x))
    else:
        return _func(x)


def get_vocab_color(x, vocab, default="black"):
    """Get a controlled vocabulary color for a value"""

    if vocab not in VOCAB_COLOR_MAPS:
        raise ValueError(f"Vocabulary `{vocab}` is not defined!")
    color = mcolors.to_rgb(VOCAB_COLOR_MAPS[vocab].get(x, default))
    return color


COLOR_MAPS = dict(
    solar_extra=[
        "#3361A5",
        "#248AF3",
        "#14B3FF",
        "#88CEEF",
        "#C1D5DC",
        "#EAD397",
        "#FDB31A",
        "#E42A2A",
        "#A31D1D",
    ],
    solar_blues=[
        "#FCFCFC",
        "#C0E4FD",
        "#75CEFE",
        "#0CB9FF",
        "#1BA7FF",
        "#1E95FF",
        "#2884E7",
        "#3072C5",
        "#3361A5",
    ],
    solar_rojos=[
        "#FCFCFC",
        "#FFEDB0",
        "#FFDF5F",
        "#FEC510",
        "#FA8E24",
        "#F14C2B",
        "#DA2828",
        "#BE2222",
        "#A31D1D",
    ],
    aqua_brick=[
        "#019bcf",
        "#0d94c2",
        "#1a8cb4",
        "#2685a7",
        "#337d99",
        "#3f768c",
        "#4c6e7e",
        "#586771",
        "#655f63",
        "#715856",
        "#7e5048",
        "#8a493b",
        "#97412d",
        "#a33a20",
    ],
    lawhoops=[
        "#371377",
        "#7700FF",
        "#9E0142",
        "#CB2314",
        "#FF0080",
        "#DC494C",
        "#F88D51",
        "#FAD510",
        "#FFFF5F",
        "#88CFA4",
        "#238B45",
        "#02401B",
        "#0AD7D3",
        "dodgerblue",
        "#046C9A",
        "#273046",
        "#A2A475",
        "#354823",
        "#595959",
        "#1E1E1E",
    ],
    # wes anderson palettes
    chevalier=["#446455", "#FDD262", "#D3DDDC", "#C7B19C"],
    zissou=["#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00"],
    fantastic_fox=["#DD8D29", "#E2D200", "#46ACC8", "#E58601", "#B40F20"],
    darjeeling=["#FF0000", "#00A08A", "#F2AD00", "#F98400", "#5BBCD6"],
    wolfgang_basic=[
        "#FFFFD9",
        "#EDF8B1",
        "#C7E9B4",
        "#7FCDBB",
        "#41B6C4",
        "#1D91C0",
        "#225EA8",
        "#253494",
        "#081D58",
    ],
    wolfgang_extra=[
        "#FFFFFF",
        "#FCFED3",
        "#E3F4B1",
        "#ABDEB6",
        "#60C1BF",
        "#2A9EC1",
        "#206AAD",
        "#243996",
        "#081D58",
    ],
)

# Register colormaps for global usage
for _name, _palette in COLOR_MAPS.items():
    color_list = list(map(mcolors.to_rgb, _palette))

    _cmap = mcolors.LinearSegmentedColormap.from_list(_name, color_list)
    _cmap_r = mcolors.LinearSegmentedColormap.from_list(_name + "_r", color_list[::-1])

    for orient, loc_cmap in zip(["", "_r"], [_cmap, _cmap_r]):
        name = _name + orient
        try:
            mcm.get_cmap(name)
        except ValueError:
            locals()[name] = loc_cmap
            try:
                mcm.register_cmap(name, loc_cmap)
            except AttributeError:
                mpl.colormaps.register(loc_cmap, name=name)

