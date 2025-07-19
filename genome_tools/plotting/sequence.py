import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from shapely.geometry import Polygon as ShapelyPolygon

from .colors.cm import VOCAB_COLOR_MAPS

from fontTools.ttLib import TTFont
from fontTools.pens.basePen import BasePen

from matplotlib import font_manager

import string


class MultiPolygonPen(BasePen):
    def __init__(self, glyphSet, approximation_scale=1):
        BasePen.__init__(self, glyphSet)
        self.polygons = []  # List to hold multiple polygons (each a list of points)
        self.current_polygon = []  # Current polygon being drawn
        self.approximation_scale = approximation_scale

    def _moveTo(self, p):
        # If moving to a new point not immediately after closing a path, start a new polygon
        if self.current_polygon:
            self.polygons.append(self.current_polygon)
            self.current_polygon = []
        self.current_polygon.append(p)

    def _lineTo(self, p):
        self.current_polygon.append(p)

    def _curveToOne(self, p1, p2, p3):
        if self.current_polygon:
            start = self.current_polygon[-1]
        else:
            start = (0, 0)

        num_segments = int(self._curve_length(start, p1, p2, p3) * self.approximation_scale)

        for i in range(1, num_segments + 1):
            t = i / num_segments
            point = self._get_bezier_point(t, start, p1, p2, p3)
            self.current_polygon.append(point)

    def _closePath(self):
        # Close the path by ensuring the last point is the same as the first
        if self.current_polygon and self.current_polygon[0] != self.current_polygon[-1]:
            self.current_polygon.append(self.current_polygon[0])
        self.polygons.append(self.current_polygon)
        self.current_polygon = []

    def _get_bezier_point(self, t, start, p1, p2, p3):
        """Calculate a point in a cubic Bezier curve."""
        x = (1 - t)**3 * start[0] + 3 * (1 - t)**2 * t * p1[0] + 3 * (1 - t) * t**2 * p2[0] + t**3 * p3[0]
        y = (1 - t)**3 * start[1] + 3 * (1 - t)**2 * t * p1[1] + 3 * (1 - t) * t**2 * p2[1] + t**3 * p3[1]
        return (x, y)

    def _curve_length(self, start, p1, p2, p3, steps=10):
        """Approximate the length of the bezier curve."""
        length = 0
        last_point = start
        for step in range(1, steps + 1):
            t = step / steps
            current_point = self._get_bezier_point(t, start, p1, p2, p3)
            length += np.linalg.norm(np.array(current_point) - np.array(last_point))
            last_point = current_point
        return length

    def get_multipolygon(self):
        # Call this method to get the multipolygon after drawing is complete
        if self.current_polygon:
            # If a polygon is in progress, add it to the list
            self._closePath()
        return self.polygons


def standardize_multipolygon(multipolygon):
    """
    Standardize the coordinates of a multipolygon so its bounding box fits in a 1x1 grid,
    maintaining the aspect ratio and centering it within the grid.

    Parameters:
    - multipolygon: A list of lists of (x, y) tuples.

    Returns:
    - A standardized multipolygon as a list of lists of (x, y) tuples.
    """
    # Flatten the list of points to calculate the bounding box
    all_points = [point for polygon in multipolygon for point in polygon]
    min_x = min(point[0] for point in all_points)
    max_x = max(point[0] for point in all_points)
    min_y = min(point[1] for point in all_points)
    max_y = max(point[1] for point in all_points)

    # Calculate scale based on the larger dimension of the bounding box
    width, height = max_x - min_x, max_y - min_y
    scale_x = 1 / width
    scale_y = 1 / height

    # Calculate translation to center the polygon in the 1x1 grid
    translate_x = -(min_x + max_x) / 2 * scale_x + 0.5
    translate_y = -(min_y + max_y) / 2 * scale_y + 0.5

    # Apply scaling and translation
    standardized_multipolygon = [
        [(x * scale_x + translate_x, y * scale_y + translate_y) for x, y in polygon]
        for polygon in multipolygon
    ]

    return standardized_multipolygon


def multipolygon_to_wkt(multipolygon):
    """
    Converts a multipolygon structure into a WKT MULTIPOLYGON or POLYGON string.

    Parameters:
    - multipolygon: A list of lists of (x, y) tuples.

    Returns:
    - A string in WKT format.
    """
    if len(multipolygon) == 1:
        # Only one polygon, no holes
        polygon_str = "POLYGON ((" + ", ".join(f"{x} {y}" for x, y in multipolygon[0]) + "))"
        return polygon_str
    else:
        # Multiple polygons, including holes
        multipolygon_str = "MULTIPOLYGON ("
        polygons_str = []
        for polygon in multipolygon:
            polygons_str.append("((" + ", ".join(f"{x} {y}" for x, y in polygon) + "))")
        multipolygon_str += ", ".join(polygons_str) + ")"
        return multipolygon_str


def get_polygon(glyph_name, glyph_set):
    pen = MultiPolygonPen(glyph_set, approximation_scale=0.03)
    glyph = glyph_set[glyph_name]
    glyph.draw(pen)
    return standardize_multipolygon(pen.get_multipolygon())


def get_letter_polygons(font_name, letters):
    if os.path.isfile(font_name):
        font = TTFont(font_name)
    else:
        font = TTFont(font_manager.findfont(font_name))
    glyph_set = font.getGlyphSet()
    return {l: get_polygon(l, glyph_set) for l in letters}

default_font = 'Arial'
default_letter_polygons = get_letter_polygons(default_font, string.ascii_uppercase + string.ascii_lowercase)


def transform_polygon(polygon, width_scale=1.0, height_scale=1.0, x_offset=0.0, y_offset=0.0):
    """
    Apply scaling to the width relative to the middle of the polygon (x=0.5) and
    apply uniform scaling to the height. Also apply translation as specified.

    Parameters:
    - polygon: A list of (x, y) tuples representing the polygon's vertices.
    - width_scale: Scaling factor in the x-direction.
    - height_scale: Scaling factor in the y-direction.
    - x_offset: Translation offset in the x-direction.
    - y_offset: Translation offset in the y-direction.

    Returns:
    - Transformed list of (x, y) tuples.
    """
    center_x = 0
    
    transformed_polygon = []
    for x, y in polygon:
        new_x = ((x - center_x) * width_scale) + center_x + x_offset
        new_y = (y * height_scale) + y_offset
        transformed_polygon.append((new_x, new_y))
    
    return transformed_polygon


def classify_polygons_by_nesting_depth(polygons):
    """
    Classify each polygon by the number of times it is contained within other polygons (its depth).
    
    Parameters:
    - polygons: A list of lists of (x, y) tuples representing the polygons.

    Returns:
    - A list of tuples, where each tuple contains (depth, polygon).
    """
    shapely_polygons = [ShapelyPolygon(p) for p in polygons]

    # List to store (depth, polygon) tuples
    ranked_polygons = []

    # Determine depth for each polygon
    for i, polygon in enumerate(shapely_polygons):
        depth = 0  # Start with depth 0

        # Check containment within all other polygons
        for j, other_polygon in enumerate(shapely_polygons):
            if i != j and polygon.within(other_polygon):
                depth += 1  # Increase depth for each containment

        ranked_polygons.append((depth, polygons[i]))

    # Sort polygons by depth
    ranked_polygons.sort(key=lambda x: x[0])

    return ranked_polygons


def add_letter_to_axis(ax, multipolygon, col, x, y, height, width_scale=1.0):
    """
    Add 'let' with position x, y and height height to matplotlib axis 'ax', adjusting width by width_scale.
    """
    # Classify polygons by their containment depth
    ranked_polygons = classify_polygons_by_nesting_depth(multipolygon)

    # Draw polygons from lowest depth to highest depth
    for depth, polygon in ranked_polygons:
        # Determine face color based on depth parity
        face_color = col if depth % 2 == 0 else 'white'
        
        # Transform and plot the polygon
        transformed_polygon = transform_polygon(polygon, width_scale, height, x, y)
        polygon_patch = Polygon(transformed_polygon, closed=True, edgecolor='none', facecolor=face_color, linewidth=0)
        ax.add_patch(polygon_patch)


def seq_plot(letter_heights, ax=None, vocab="dna", offset=0, width_scale=1.0, font=default_font, **kwargs):
    if font == default_font:
        letter_polygons = default_letter_polygons
    else:
        letter_polygons = get_letter_polygons(font, string.ascii_uppercase + string.ascii_lowercase)

    if not ax:
        ax = plt.gca()
    fig = ax.figure

    assert letter_heights.shape[1] == len(VOCAB_COLOR_MAPS[vocab])
    x_range = [0, letter_heights.shape[0]]
    pos_heights = np.copy(letter_heights)
    pos_heights[letter_heights < 0] = 0
    neg_heights = np.copy(letter_heights)
    neg_heights[letter_heights > 0] = 0
    
    max_pos_h = 0.0
    min_neg_h = 0.0

    for x_pos, heights in enumerate(letter_heights):
        letters_and_heights = sorted(zip(heights, list(VOCAB_COLOR_MAPS[vocab].keys())))
        y_pos_pos = 0.0
        y_neg_pos = 0.0
        # x_pos += interval.start
        for height, letter in letters_and_heights:
            color = VOCAB_COLOR_MAPS[vocab][letter]
            polygons = letter_polygons[letter]
            if height > 0:
                add_letter_to_axis(
                    ax, polygons, color, x_pos + offset, y_pos_pos, height,
                    width_scale=width_scale
                )
                y_pos_pos += height
            elif height < 0:
                add_letter_to_axis(
                    ax, polygons, color, x_pos + offset, y_neg_pos, height,
                    width_scale=width_scale
                )
                y_neg_pos += height
                
        max_pos_h = max(max_pos_h, y_pos_pos)
        min_neg_h = min(min_neg_h, y_neg_pos)

    ax.set_xlim(left=offset, right=len(letter_heights) + offset)
    ax.set_ylim(bottom=min_neg_h, top=max_pos_h)
