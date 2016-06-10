# Copyright 2016 Jeff Vierstra

import sys

from matplotlib.ticker import MaxNLocator, FuncFormatter

class track(object):
    """Base level track object"""
    def __init__(self, interval, **kwargs):
    
        self.interval = interval

        self.options = {
            'edge_color': 'none',
            'lw': 1
        }

        self.options.update(kwargs)

    def format_axis(self, ax):
        
        [spine.set_color('none') for loc, spine in ax.spines.items()]
        ax.xaxis.set_tick_params(direction = 'out')
        ax.yaxis.set_tick_params(direction = 'out')
        
        # set and format xlim
        ax.set_xlim(left = 0, right = len(self.interval))
        
        formatter = FuncFormatter(lambda x, pos: int(x + self.interval.start))
        ax.xaxis.set(major_formatter = formatter)
        
        ax.xaxis.set(visible = False)
        ax.yaxis.set(visible = False)

    def render(self, ax):
        
         pass

class row_element(object):

    def __init__(self, interval, prev = None, next = None):
        self.interval = interval
        self.start = interval.start - 5
        self.end = interval.end + 5
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

def pack_rows(intervals):
    
    rows = []
    row_num = {}
    
    nrows = 0

    for interval in intervals:
        #sys.stderr.write("adding feat...\n")
        re = row_element(interval)
        
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
                raise Exception("could not place element in new row")
                
    return nrows, row_num

