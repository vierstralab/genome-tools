import os, string
import numpy as np

import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.path import Path
from matplotlib.patches import PathPatch

from fontTools.ttLib import TTFont
from fontTools.pens.basePen import BasePen

from shapely.geometry import Polygon as ShapelyPolygon, MultiPolygon, GeometryCollection
from shapely.ops import unary_union
from shapely.geometry.polygon import orient

from genome_tools.plotting.colors.cm import VOCAB_COLOR_MAPS
from genome_tools.data.pwm import relative_info_content


class RingPen(BasePen):
    """
    Records glyph outlines as a list of rings (each ring = list of (x,y) points, closed).
    Cubic curves are approximated by line segments (no forced min segment count).
    """
    def __init__(self, glyphSet, approximation_scale=0.03):
        super().__init__(glyphSet)
        self.rings = []
        self.cur = []
        self.approximation_scale = float(approximation_scale)

    def _moveTo(self, p):
        if self.cur:
            self.rings.append(self.cur)
            self.cur = []
        self.cur = [p]

    def _lineTo(self, p):
        self.cur.append(p)

    def _get_bezier_point(self, t, start, p1, p2, p3):
        x = (1 - t)**3 * start[0] + 3 * (1 - t)**2 * t * p1[0] + 3 * (1 - t) * t**2 * p2[0] + t**3 * p3[0]
        y = (1 - t)**3 * start[1] + 3 * (1 - t)**2 * t * p1[1] + 3 * (1 - t) * t**2 * p2[1] + t**3 * p3[1]
        return (x, y)

    def _curve_length(self, start, p1, p2, p3, steps=10):
        length = 0.0
        last = start
        for k in range(1, steps + 1):
            t = k / steps
            cur = self._get_bezier_point(t, start, p1, p2, p3)
            length += np.linalg.norm(np.array(cur) - np.array(last))
            last = cur
        return length

    def _curveToOne(self, p1, p2, p3):
        start = self.cur[-1] if self.cur else (0.0, 0.0)
        num_segments = int(self._curve_length(start, p1, p2, p3) * self.approximation_scale)

        # IMPORTANT: per your request, no max(1, ...) here.
        # If num_segments is 0, just add the endpoint once (prevents div-by-zero).
        if num_segments <= 0:
            self.cur.append(p3)
            return

        for i in range(1, num_segments + 1):
            t = i / num_segments
            self.cur.append(self._get_bezier_point(t, start, p1, p2, p3))

    def _closePath(self):
        if self.cur and self.cur[0] != self.cur[-1]:
            self.cur.append(self.cur[0])
        if self.cur:
            self.rings.append(self.cur)
        self.cur = []

    def get_rings(self):
        if self.cur:
            self._closePath()
        return self.rings


def standardize_rings(rings):
    """
    Uniform scale + center so bbox fits in [0,1]^2 preserving aspect ratio.
    """
    pts = [p for r in rings for p in r]
    if not pts:
        return rings
    xs, ys = zip(*pts)
    minx, maxx = min(xs), max(xs)
    miny, maxy = min(ys), max(ys)
    w, h = maxx - minx, maxy - miny
    if w == 0 or h == 0:
        return rings

    s = 1.0 / max(w, h)
    cx, cy = (minx + maxx) / 2.0, (miny + maxy) / 2.0
    tx, ty = 0.5 - s * cx, 0.5 - s * cy
    return [[(s * x + tx, s * y + ty) for x, y in r] for r in rings]


def signed_area(ring):
    """
    Signed area: >0 means CCW, <0 means CW (for typical coordinate systems).
    Works if ring is closed or not.
    """
    n = len(ring)
    if n < 3:
        return 0.0
    s = 0.0
    for i in range(n - 1):
        x1, y1 = ring[i]
        x2, y2 = ring[i + 1]
        s += x1 * y2 - x2 * y1
    if ring[0] != ring[-1]:
        x1, y1 = ring[-1]
        x2, y2 = ring[0]
        s += x1 * y2 - x2 * y1
    return 0.5 * s


def rings_to_geometry_winding(rings, repair_invalid=True, try_swap_if_empty=True):
    """
    Build a SINGLE geometry from rings using winding:
        G = union(pos_winding) - union(neg_winding)

    Returns Polygon/MultiPolygon/GeometryCollection.
    """
    pos, neg = [], []
    for r in rings:
        p = ShapelyPolygon(r)
        if repair_invalid and (not p.is_valid):
            p = p.buffer(0)
        if p.is_empty:
            continue
        (pos if signed_area(r) > 0 else neg).append(p)

    if not pos and not neg:
        return GeometryCollection()

    Upos = unary_union(pos) if pos else GeometryCollection()
    Uneg = unary_union(neg) if neg else GeometryCollection()

    g = Upos.difference(Uneg)

    # Some fonts flip winding convention; optionally try swapping once.
    if try_swap_if_empty and g.is_empty and (not Upos.is_empty or not Uneg.is_empty):
        g2 = Uneg.difference(Upos)
        if (not g2.is_empty) and (g2.area > g.area):
            g = g2

    # Normalize for stable export/path building (outer CCW, holes CW)
    if g.geom_type == "Polygon":
        g = orient(g, sign=1.0)
    elif g.geom_type == "MultiPolygon":
        g = MultiPolygon([orient(p, sign=1.0) for p in g.geoms])

    return g


def get_glyph_geometry(char, glyph_set, approximation_scale=0.03):
    pen = RingPen(glyph_set, approximation_scale=approximation_scale)
    glyph_set[char].draw(pen)
    rings = standardize_rings(pen.get_rings())
    return rings_to_geometry_winding(rings)


def get_letter_geometries(font_name, letters, approximation_scale=0.03):
    """
    dict: char -> shapely geometry
    """
    font = TTFont(font_name) if os.path.isfile(font_name) else TTFont(font_manager.findfont(font_name))
    glyph_set = font.getGlyphSet()
    out = {}
    for ch in letters:
        if ch in glyph_set:
            out[ch] = get_glyph_geometry(ch, glyph_set, approximation_scale=approximation_scale)
    return out


default_font = "Arial"
default_letter_geoms = get_letter_geometries(default_font, string.ascii_uppercase + string.ascii_lowercase)


def _ring_to_path(ring):
    pts = list(ring)
    if not pts:
        return [], []
    if pts[0] != pts[-1]:
        pts.append(pts[0])

    verts = [(pts[0][0], pts[0][1])]
    codes = [Path.MOVETO]
    for x, y in pts[1:]:
        verts.append((x, y))
        codes.append(Path.LINETO)

    # CLOSEPOLY wants a final vertex
    verts.append((pts[0][0], pts[0][1]))
    codes.append(Path.CLOSEPOLY)
    return verts, codes


def geometry_to_path(geom):
    """
    Convert Shapely Polygon/MultiPolygon into one compound Path:
    includes exteriors and interiors (holes).
    """
    if geom.is_empty:
        return Path([], [])

    verts, codes = [], []

    def add_polygon(poly):
        nonlocal verts, codes
        v, c = _ring_to_path(list(poly.exterior.coords))
        verts += v; codes += c
        for interior in poly.interiors:
            v, c = _ring_to_path(list(interior.coords))
            verts += v; codes += c

    if geom.geom_type == "Polygon":
        add_polygon(geom)
    elif geom.geom_type == "MultiPolygon":
        for poly in geom.geoms:
            add_polygon(poly)
    else:
        return Path([], [])

    return Path(verts, codes)


def transform_path(path, width_scale=1.0, height_scale=1.0, x_offset=0.0, y_offset=0.0):
    """
    Standardized glyph is in ~[0,1]x[0,1].
    Scale x about 0.5, scale y about 0 (same behavior as your original y*height + y_offset).
    """
    if path.vertices.size == 0:
        return path
    v = path.vertices.copy()
    v[:, 0] = (v[:, 0] - 0.5) * width_scale + 0.5 + x_offset
    v[:, 1] = v[:, 1] * height_scale + y_offset
    return Path(v, path.codes)


def add_geometry_to_axis(ax, geom, color, x, y, height, width_scale=1.0, center_scale=True):
    if center_scale:
        x_off = x + (1 - width_scale) / 2
    else:
        x_off = x

    path = geometry_to_path(geom)
    path = transform_path(path, width_scale=width_scale, height_scale=height, x_offset=x_off, y_offset=y)

    patch = PathPatch(path, facecolor=color, edgecolor="none", linewidth=0)

    ax.add_patch(patch)
    return patch


def plot_letter(letter, x, y, height=1.0, width=1.0, center_scale=True,
                vocab="dna", color=None, font=default_font, ax=None):
    if ax is None:
        ax = plt.gca()

    geoms = default_letter_geoms if font == default_font else get_letter_geometries(
        font, string.ascii_uppercase + string.ascii_lowercase
    )

    geom = geoms[letter]
    if color is None:
        color = VOCAB_COLOR_MAPS[vocab].get(letter, "black")

    add_geometry_to_axis(ax, geom, color, x, y, height, width_scale=width, center_scale=center_scale)
    return ax


def seq_plot(letter_heights: np.ndarray, ax=None, vocab="dna", offset=0,
             width_scale=1.0, font=default_font, center_scale=True):
    geoms = default_letter_geoms if font == default_font else get_letter_geometries(
        font, string.ascii_uppercase + string.ascii_lowercase
    )

    if ax is None:
        ax = plt.gca()

    if isinstance(vocab, str):
        vocab = VOCAB_COLOR_MAPS[vocab]

    letters = list(vocab.keys())
    assert letter_heights.shape[1] == len(letters)

    max_pos_h = 0.0
    min_neg_h = 0.0

    for x_pos, heights in enumerate(letter_heights):
        y_pos = 0.0
        y_neg = 0.0

        for h, letter in sorted(zip(heights, letters)):
            if h == 0:
                continue
            geom = geoms[letter]
            col = vocab[letter]

            if h > 0:
                add_geometry_to_axis(ax, geom, col, x_pos + offset, y_pos, h,
                                     width_scale=width_scale, center_scale=center_scale)
                y_pos += h
            else:
                add_geometry_to_axis(ax, geom, col, x_pos + offset, y_neg, h,
                                     width_scale=width_scale, center_scale=center_scale)
                y_neg += h

        max_pos_h = max(max_pos_h, y_pos)
        min_neg_h = min(min_neg_h, y_neg)

    ax.set_xlim(left=offset, right=len(letter_heights) + offset)
    ax.set_ylim(bottom=min_neg_h, top=max_pos_h)
    return ax


def plot_motif_logo(pfm: np.ndarray, ax=None, offset=-0.5, rc=False, **kwargs):
    if rc:
        pfm = pfm[::-1, ::-1]
    seq_plot(relative_info_content(pfm.T), ax=ax, offset=offset, **kwargs)
    return ax
