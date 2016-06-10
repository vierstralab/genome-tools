# Copyright 2016 Jeff Vierstra

import numpy as np

from matplotlib.ticker import MaxNLocator

from track import track

class continuous_data_track(track):
    
    def __init__(self, interval, data, **kwargs):
        super(continuous_data_track, self).__init__(interval, **kwargs)
     
        self.data = data

    def format_axis(self, ax):
        super(continuous_data_track, self).format_axis(ax)

        ax.spines['bottom'].set_color('black')

        locator = MaxNLocator(4, prune = 'upper')
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

    def render(self, ax):
    
        self.format_axis(ax)
    
        xx = self.step(np.arange(len(self.interval)), xaxis = True)
        yy = self.step(self.data)
    
        if 'fill_between' in self.options:
            ax.fill_between(xx, np.ones(len(xx)) * self.options['fill_between'], yy, facecolor = self.options.get('fill_color', 'blue'), edgecolor = self.options['edge_color'])
        else:
            ax.plot(xx, yy)
        
        #self.set_xaxis(ax)
        