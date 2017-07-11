# Copyright 2016 Jeff Vierstra

import numpy as np

from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import blended_transform_factory

from track import track, segment

class continuous_data_track(track):
    
    def __init__(self, interval, data, **kwargs):
        super(continuous_data_track, self).__init__(interval, **kwargs)
     
        self.data = data

    def format_axis(self, ax):
        super(continuous_data_track, self).format_axis(ax)

        ax.spines['bottom'].set_color('black')

        locator = MaxNLocator(3, prune = 'upper')
        ax.yaxis.set(visible = True, major_locator = locator)

        if 'min' in self.options:
            ax.set_ylim(bottom = self.options['min'])

        if 'max' in self.options:
            ax.set_ylim(top = self.options['max'])

    def step(self, vals, xaxis = False, step_interval = 0):
        if xaxis and step_interval == 0:
            step_interval = abs(vals[1] - vals[0]) / 2.0
        step_vals = np.array(zip(vals - step_interval, vals + step_interval)).ravel()
        return step_vals

    def simplify(self, step = 100):
        xx = []
        yy = []
        for i in range(step, len(self.data)-step, step):
            xx.append(self.interval.start + i)
            yy.append(np.amax(self.data[i-step:i+step]))
        return np.array(xx), np.array(yy)

    def render(self, ax):
    
        self.format_axis(ax)

        if 'simplify' in self.options:
            xx, yy = self.simplify(step = self.options['simplify'])
            xs = self.step(xx, xaxis = True, step_interval = 100)
            ys = self.step(yy, step_interval = 0)
        else:
            xs = self.step(np.arange(self.interval.start, self.interval.end), xaxis = True)
            ys = self.step(self.data)
    
        (ybot, ytop) = ax.get_ylim()
        ys[ys > ytop] = ytop

        if 'fill_between' in self.options:
            ax.fill_between(xs, np.ones(len(xs)) * self.options['fill_between'], ys, facecolor = self.options.get('fill_color', 'blue'), edgecolor = self.options['edge_color'])
        else:
            ax.plot(xs, ys)

        if 'clip_bar' in self.options:
            trans = blended_transform_factory(ax.transData, ax.transAxes)
            for x0, x1 in segment(self.data, ytop):
                ax.plot([x0 + self.interval.start, x1 + self.interval.start], [1, 1], color = self.options['clip_bar']['color'], lw = self.options['clip_bar']['lw'], transform = trans, clip_on = False)

        