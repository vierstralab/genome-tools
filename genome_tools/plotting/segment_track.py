import numpy as numpy
import pysam

from .. import bed
from .. import genomic_interval

from matplotlib.patches import Rectangle
from matplotlib.transforms import Affine2D

from track import track, pack_rows


class segment_track(track):

	def __init__(self, interval):
		super(segment_track, self).__init__(interval)

		self.segments = []

	def load_from_file(self, file_path, file_type = 'tabix', file_format = 'bed3'):

		tabix = pysam.TabixFile(file_path)
		data = tabix.fetch(self.interval.chrom, self.interval.start, self.interval.end, parser = None)

		data_iterator = bed.bed3_iterator(data) if file_format == 'bed3' else bed.bed6_iterator(data)

		for segment in data_iterator:
			self.segments.append(segment)

	def format_axis(self, ax):

		super(segment_track, self).format_axis(ax)

		ax.xaxis.set(visible = False)
		ax.yaxis.set(visible = False)

	def render(self, ax, pack = True, labels = False):

		self.format_axis(ax)

		nrows, rows = pack_rows(self.segments)
		
		for segment in self.segments:

			n = rows[segment] if pack else 1
			
			if self.orientation == 'horizontal':
				
				x0 = segment.start
				x1 = segment.end
				y0 = -0.2+n
				w = x1-x0
				h = 0.4

				label_x, label_y, label_rot = x1, y0, 0

			else:
				
				y0 = segment.start
				y1 = segment.end
				x0 = -0.2+n
				w = 0.4
				h = y1-y0

				label_x, label_y, label_rot = x0, y1, 90

			p = Rectangle((x0, y0), w, h, edgecolor = 'none', facecolor = 'goldenrod', zorder = 1)
			ax.add_patch(p)
			
			if labels:
				ax.text(label_x, label_y, segment.name, rotation = label_rot, ha = 'left', va = 'bottom')

		#ax.set_ylim(bottom = 0.5, top = (nrows if pack else 1) + 1.5)
		if self.orientation == 'horizontal':
			ax.set_ylim(bottom = 0.5, top = (nrows if pack else 1) + 1.5)
		else:
			ax.set_xlim(left = 0.5, right = (nrows if pack else 1) + 1.5)

			
