# Copyright 2016 Jeff Vierstra

import numpy as np

from track import track

class heatmap_data_track(track):

	def __init__(self, interval, data, **kwargs):
		super(heatmap_data_track, self).__init__(interval, **kwargs)
     
		self.data = data

	def format_axis(self, ax):
		super(heatmap_data_track, self).format_axis(ax)

		ax.yaxis.set(visible = True, ticks = np.arange(self.data.shape[0]), ticklabels = [])
		#ax.set_ylim([0, self.data.shape[0]])

	def render(self, ax):
		
		self.format_axis(ax)

		ax.imshow(self.data, aspect = "auto", interpolation = "nearest")
