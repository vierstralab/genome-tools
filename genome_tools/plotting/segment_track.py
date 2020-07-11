import numpy as numpy
import pysam


from .. import load_segmented_data

from matplotlib.patches import Rectangle
from matplotlib.transforms import Affine2D

from .track import track, pack_rows


class segment_track(track):

	def __init__(self, interval, **kwargs):
		super(segment_track, self).__init__(interval, **kwargs)

		self.segments = []

	def load_data(self, filepath, filterfunc=None):
		segments = load_segmented_data(filepath, self.interval)
		self.segments = list(filter(filterfunc, segments)) if filterfunc else segments

	def format_axis(self, ax):
		super(segment_track, self).format_axis(ax)

		ax.xaxis.set(visible = False)
		ax.yaxis.set(visible = False)

	def render(self, ax, pack = True, labels = False):

		self.format_axis(ax)
		self.format_spines(ax, remove_spines=['top', 'right'])

		nrows, rows = pack_rows(self.segments)
		
		for segment in self.segments:

			n = rows[segment] if pack else 1
							
			x0 = segment.start
			x1 = segment.end
			y0 = -0.2+n
			w = x1-x0
			h = 0.4

			label_x, label_y, label_rot = x1, y0, 0

			p = Rectangle((x0, y0), w, h, 
				edgecolor=self.options["edgecolor"], 
				facecolor=self.options["facecolor"],
				lw=self.options["lw"],
				ls=self.options["ls"],
				zorder=1)

			ax.add_patch(p)
			
			if labels:
				ax.text(label_x, label_y, segment.name, rotation=label_rot, ha='left', va='bottom')

		ax.set_ylim(bottom=0.5, top=(nrows if pack else 1) + 1.5)
		ax.invert_yaxis()

		super(segment_track, self).render(ax)
