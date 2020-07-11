# Copyright 2019 Jeff Vierstra

# Draw an ideogram; one can download the ideogram files from the UCSC genome browser

from ..helpers import open_file


from matplotlib.collections import BrokenBarHCollection
from matplotlib.patches import Polygon

import numpy as np

color_lookup = {
                  'gneg': (1., 1., 1.),
                'gpos25': (.6, .6, .6),
                'gpos50': (.4, .4, .4),
                'gpos75': (.2, .2, .2),
               'gpos100': (0., 0., 0.),
                  'acen': (.8, .4, .4),
                  'gvar': (.8, .8, .8),
                 'stalk': (.9, .9, .9),
               }

class ideogram(object):

    def __init__(self, filepath):
        self.xranges = None
        self.colors = None
        self.centromeres = None

        self.read_ideogram(filepath)

    def read_ideogram(self, filepath):
        
        xranges={}
        colors={}
        centromeres={}

        f = open_file(filepath)

        last_chrom=None
        xr=[]
        for line in f:
            chrom, start, stop, label, stain = line.strip().split('\t')
            start = int(start)
            stop = int(stop)
            width = stop - start
        
            if stain=="acen":
                centromeres[chrom] = centromeres.get(chrom, []) + [(start, width)]
                continue
        
            xranges[chrom] = xranges.get(chrom, []) + [(start, width)]
            colors[chrom] = colors.get(chrom, []) + [(color_lookup[stain])]

        f.close()

        self.xranges = xranges
        self.colors = colors
        self.centromeres = centromeres

    def render(self, ax, chrom, pos = None):
        
        try:
            xranges = self.xranges[chrom]
            colors = self.colors[chrom]
            centromeres = self.centromeres[chrom]
        except:
            print("Error: No chromosome named: {}".format(chrom))
            return

        yranges = (0, 0.5)
        
        coll = BrokenBarHCollection(xranges, yranges, facecolors=colors, edgecolors='black', linewidths=0.5)    
        ax.add_collection(coll)

        if pos:
            ax.axvline(pos, color='red', lw=4)
        w = xranges[-1][0]+xranges[-1][1]
        
        pad=w*0.05
        
        ax.set_xlim(0-pad, xranges[-1][0]+xranges[-1][1]+pad)
        ax.xaxis.set_visible(False)
        
        center = yranges[0] + yranges[1]/2.
        
     
        x0, y0 = centromeres[0][0], yranges[0]
        x1, y1 = centromeres[0][0], yranges[1]
        x2, y2 = centromeres[0][0]+centromeres[0][1], center
        cent = Polygon(np.array([[x0, y0], [x1, y1], [x2, y2]]), closed=True, 
                       fc=color_lookup['acen'], ec='black', linewidth=0.5)
        ax.add_patch(cent)
        
        x0, y0 = centromeres[1][0], center
        x1, y1 = centromeres[1][0]+centromeres[1][1], yranges[1]
        x2, y2 = centromeres[1][0]+centromeres[1][1], yranges[0]
        
        cent = Polygon(np.array([[x0, y0], [x1, y1], [x2, y2]]), closed=True,
                        fc=color_lookup['acen'], ec='black', linewidth=0.5)
        ax.add_patch(cent)
        
        
        ax.set_yticks([center])
        ax.set_yticklabels([chrom])
        ax.set_ylim(-0.2, 0.7)
        
        [ax.spines[loc].set_color('none') for loc in ['top', 'left','right','bottom']]