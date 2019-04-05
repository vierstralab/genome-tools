# Copyright 2016 Jeff Vierstra

import numpy as np

from .. import load_data
from .track import track

class heatmap_track(track):

	def __init__(self, interval, data, **kwargs):
		super(heatmap_track, self).__init__(interval, **kwargs)

		self.data = data if data.ndim>1 else data[np.newaxis,:]

	def load_data(self, filepath, column=5, dtype=float):
		self.data = load_data(filepath, interval, [column], dtype).T

	def format_axis(self, ax):
		super(heatmap_track, self).format_axis(ax)

		ax.yaxis.set(visible = True, ticks = np.arange(0.5, self.data.shape[0], 0.5), ticklabels = [])
		#ax.set_ylim([0, self.data.shape[0]])

	def render(self, ax):
		
		self.format_axis(ax)
		self.format_spines(ax, remove_spines=[])

		x = range(self.interval.start, self.interval.end+1)
		y = range(self.data.shape[0]+1)
		
		vmin = self.options.get('vmin', np.min(self.data))
		vmax = self.options.get('vmax', np.max(self.data))
		vmid = self.options.get('vmid', 0)

		if vmax>0 and vmin<0:
			divergent_data=True
			vmin+=vmid
			vmax+=vmid

		if not 'cmap' in self.options:
			if divergent_data:
				cmap='RdBu'
			else:
				cmap='RdPu'
		else:
			cmap = self.options['cmap']

		ax.pcolormesh(x, y, self.data[::-1,:], cmap = cmap, vmin=vmin, vmax=vmax, rasterized=True)

		super(heatmap_track, self).render(ax)