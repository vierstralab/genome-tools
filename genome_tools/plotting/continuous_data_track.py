# Copyright 2016 Jeff Vierstra

import numpy as np

from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import blended_transform_factory

from .. import load_data
from .track import track, segment

class continuous_data_track(track):

	def __init__(self, interval, data=None, **kwargs):
		super(continuous_data_track, self).__init__(interval, **kwargs)
		self.data = data

	def load_data(self, filepath, column=5, dtype=float):
		self.data = load_data(filepath, interval, [column], dtype)[:,0]


	def density(self, vals, window_size = 150, step_size = 20):
		"""Smooth data in density windows"""
		
		pos = np.arange(0, len(vals)-window_size, step_size)
		s = np.zeros(len(pos))
		for i, j in enumerate(pos):
			s[i] = np.sum(vals[j:j+window_size])
		return pos + (window_size/2) + self.interval.start, s

	def step(self, vals, xaxis = False, step_interval = 0):
		"""Creates a step-style plot"""

		if xaxis and step_interval == 0:
			step_interval = abs(vals[1] - vals[0]) / 2.0
		step_vals = np.array(list(zip(vals - step_interval, vals + step_interval))).ravel()
		return step_vals

	def format_axis(self, ax):
		super(continuous_data_track, self).format_axis(ax)

		locator = MaxNLocator(3, prune = 'upper')
		ax.yaxis.set(visible = True, major_locator = locator)

		if 'min' in self.options:
			ax.set_ylim(bottom = self.options['min'])

		if 'max' in self.options:
			ax.set_ylim(top = self.options['max'])


	def render(self, ax):

		if self.data is None:
			raise Exception("No data loaded!")

		self.format_axis(ax)
		self.format_spines(ax, remove_spines=['top', 'right'])

		if 'density' in self.options:
			xx, yy = self.density(self.data, window_size = self.options['density']['window_size'])
			xs = self.step(xx, xaxis = True)
			ys = self.step(yy) / self.options['density']['norm_factor']
		else:
			xs = self.step(np.arange(self.interval.start, self.interval.end), xaxis = True)
			ys = self.step(self.data)

		(ybot, ytop) = ax.get_ylim()
		ys[ys > ytop] = ytop
		ys[ys < ybot] = ybot

		if 'fill_between' in self.options:
			ax.fill_between(xs, np.ones(len(xs)) * self.options["fill_between"], ys, 
				facecolor=self.options["facecolor"], 
				edgecolor=self.options["edgecolor"],
				linewidth=self.options["lw"],
				linestyle=self.options["ls"])
		else:
			ax.plot(xs, y, 
				color=self.options["edgecolor"], 				
				lw=self.options["lw"],
				ls=self.options["ls"])

		# change to linecollection
		#if 'clip_bar' in self.options:
		#	trans = blended_transform_factory(ax.transData, ax.transAxes)
		#	for x0, x1 in segment(self.data, ytop):
		#		ax.plot([x0 + self.interval.start, x1 + self.interval.start], [1, 1], color = self.options['clip_bar']['color'], lw = self.options['clip_bar']['lw'], transform = trans, clip_on = False)

		super(continuous_data_track, self).render(ax)

