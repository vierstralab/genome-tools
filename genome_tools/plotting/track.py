# Copyright 2016 Jeff Vierstra

import sys

from matplotlib.ticker import MaxNLocator, FuncFormatter
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

class track(object):
    """Base level track object"""
    def __init__(self, interval, **kwargs):
        self.interval = interval
        self.orientation = 'horizontal'
        self.options = {
            'facecolor': 'b',
            'edgecolor': 'none',
            'lw': 1,
            'ls': '-',
            'scale_bar': None,
            'scale_bar_loc': 1
        }

        self.options.update(kwargs)

    def format_axis(self, ax):

        # set defaults
        ax.ticklabel_format(style='plain', axis = 'both', useOffset=False)
                
        # tick direction out
        ax.xaxis.set_tick_params(direction = 'out')
        ax.yaxis.set_tick_params(direction = 'out')
        
        # set x limits to interval
        ax.set_xlim(left=self.interval.start, right=self.interval.end)
        
        # x-axis choose ~3 tick positions
        locator = MaxNLocator(3, prune = 'both')
        ax.xaxis.set(major_locator = locator)
        
        # hide both x and y axis
        ax.xaxis.set(visible = False)
        ax.yaxis.set(visible = False)


    def format_spines(self, ax, remove_spines):
        all_spines = ['top', 'bottom', 'left', 'right']
        for spine in remove_spines:
            try:
                ax.spines[spine].set_visible(False)
            except KeyError:
                pass

        for spine in set(all_spines).difference(set(remove_spines)):
            try:
                ax.spines[spine].set_linewidth(0.5)
            except:
                pass

    def render(self, ax):
        # Remove the white patch behind each axes
        ax.patch.set_facecolor('none')

        # Add scale bar -- code is placed here so a scale bar can be drawn on any track instance
        if self.options['scale_bar'] is not None:
            bar = AnchoredSizeBar(ax.transData,
                self.options['scale_bar'], 
                label="%d nt" % self.options['scale_bar'], 
                loc=self.options['scale_bar_loc'],
                frameon=False)
            ax.add_artist(bar)

    def load_data(self, filepath):
        pass

class row_element(object):

    def __init__(self, interval, prev = None, next = None, padding = 5):
        self.interval = interval
        self.start = interval.start - padding
        self.end = interval.end + padding
        self.prev = prev
        self.next = next
        
class row(object):
    """ """
    def __init__(self, num):
        self.num = num
        self.elements = []
        self.first = None
        self.last = None
        
    def add_element(self, elem):
        
        if self.first is None:
            #sys.stderr.write("added feat %d-%d as only element of row %d\n" %
            #                (elem.start, elem.end, self.num))
            elem.prev = None
            elem.next = None
            self.first = self.last = elem
            return True
        
        cur = self.first
        while (cur and cur.start < elem.start):
            cur = cur.next
            
        if cur is None:
            if self.last.end < elem.start:
                #sys.stderr.write("added feat %d-%d to end of row %d\n" %
                #                 (elem.start, elem.end, self.num))
                elem.prev = self.last
                elem.next = None
                self.last.next = elem
                self.last = elem
                return True
            else:
                return False
        
        prev = cur.prev
        if prev is None:
            if elem.end < cur.start:
                #sys.stderr.write("added feat %d-%d to start of row %d "
                #                 "(before %d-%d)\n" %
                #                 (elem.start, elem.end, self.num,
                #                 cur.start, cur.end))
                elem.prev = None
                elem.next = cur
                cur.prev = elem
                self.first = elem
                return True
            else:
                return False
        
        if prev.end < elem.start and cur.start > elem.end:
            #sys.stderr.write("added feat %d-%d between %d-%d and %d-%d "
            #                 "in row %d\n" %
            #                 (elem.start, elem.end,
            #                  prev.start, prev.end,
            #                  cur.start, cur.end, self.num))
            elem.prev = prev
            elem.next = cur
            prev.next = elem
            cur.prev = elem
            return True
        
        return False

def pack_rows(intervals, padding = 5):
    
    rows = []
    row_num = {}
    
    nrows = 0

    for interval in intervals:
        #sys.stderr.write("adding feat...\n")
        re = row_element(interval, padding = padding)
        
        placed = False
        
        for r in rows:
            if(r.add_element(re)):
                placed = True
                row_num[interval] = r.num
                break
        
        if not placed:
            nrows = nrows + 1
            r = row(num = nrows)
            rows.append(r)
            row_num[interval] = r.num
            if not r.add_element(re):
                raise Exception("Could not place element in new row!")
                
    return nrows, row_num

def segment(x, threshold, w = 0, decreasing = False):
    """Segment an array into continuous elements passing a threshhold

    :returns: [(int, int), ...]: start and end points to regions that pass a threshold
    """
    dir = -1 if decreasing else 1

    ret = []
    curr_start = -1
    
    for i in range(x.shape[0]):
        if curr_start < 0:
            if dir*x[i] >= dir*threshold:
                curr_start = i-w+1
        else:
            if dir*x[i] < dir*threshold:
                if len(ret) > 0 and curr_start <= ret[-1][1]:
                    ret[-1][1] = i-1+w
                else:
                    ret.append( [curr_start, i-1+w] )
                curr_start = -1
    return ret
