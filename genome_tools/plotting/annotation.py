import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.transforms as mtransforms

def argsort(a):
    return [x for x,y in sorted(enumerate(a), key=lambda x: x[1][0])]

def bbox_expand(bbox, padding=5):
    pts = bbox.get_points()
    pts[0,:] -= padding
    pts[1,:] += padding
    return mtransforms.Bbox(pts)

def annotate(ax, labels, xy, xytext=(0, 25), slide=None, padding=5, arrowprops=None, **kwargs):
    #slide should be relevant edge of bbox - e.g. (0,0) for left, (0,1) for bottom ...
    #slide = kwargs.pop('slide', None)

    pixel_diff = 1
  
    if slide == 'h' or slide == None:
        slide0 = (0,0) # slide right
        slide1 = (1,0) # slide left
    elif slide=='v':
        slide0 = (0,1) # slide up
        slide1 = (1,1) # slide down
    else:
        raise ValueError
    
    kwargs.update(
        {
            'ha': 'center',
            'va': 'bottom'
        }
    )

    if xytext[1] < 0:
        kwargs.update({'va': 'top'})

    # Sort annotations by x axis
    sorted_idxs = argsort(xy)
    xy = [xy[i] for i in sorted_idxs]
    labels = [labels[i] for i in sorted_idxs]

    # Add annotations and draw canvas to find out where 
    # they will land on final fig
    artists = []            
    for label, data_xy in zip(labels, xy):
        a = ax.annotate(label, xy=data_xy,
            textcoords='offset pixels', 
            xytext=xytext, **kwargs)
        artists.append(a)
    plt.gcf().canvas.draw()
   
    # Split the annotations in half
    mid = len(artists)//2
    order = list(range(mid)[::-1]) + list(range(mid, len(artists)))

    #
    bboxes = []

    # Enumerate each annotation and detect collisions
    # move untill no overlaps
    for i in order:
        t = artists[i]

        # Expand patch bbox to reflect required padding
        if t.get_bbox_patch() is not None:
            curr_patch = t.get_bbox_patch()
        else:
            curr_patch = t

        curr_bbox = bbox_expand(curr_patch.get_window_extent(t._renderer), padding)

        if i >= len(artists)//2:
            slide = slide1
        else:
            slide = slide0

        # direction to slide; 
        # (0,0) = -1
        # (0,1) = -1
        # (1,0) = 1
        # (1,1) = 1
        dir  = int((slide[0] - 0.5)*2)

        # Current position; initialzed as opposite of dir
        # slide = (0, 0)
        # dir = -1
        # current = -1 * -1 * inf = inf
        current = -dir * float("inf")
    
        while 1 and slide:

            overlaps = False

            for box in bboxes:
                # .get_points() returns 2x2 array; slide selects with 'corner' to use
                box_pt = box.get_points()[slide]
                if curr_bbox.overlaps(box):
                    # dir * box_pt = -1 * x0
                    # dir * current = -1 * inf
                    # if -x0 > -inf
                    if dir * box_pt > dir * current:
                        overlaps = True
                        current = box_pt # current box x0 + padding
                        # set shift
                        # 1-slide[0], slide[1] = 1, 0 = x1
                        shift = dir * (current - curr_bbox.get_points()[1-slide[0],slide[1]])
                
            if not overlaps: 
                break

            pos = np.array(t.get_position())
            pos[slide[1]] += shift * dir * pixel_diff
            t.set_position(pos)
            
            # Redraw canvas
            plt.gcf().canvas.draw()

            # Update box; expanded to include a bit of padding
            curr_bbox = bbox_expand(curr_patch.get_window_extent(t._renderer), padding)

        bboxes.append(curr_bbox)

        if arrowprops is not None:
            pt_data = (xy[i][0], xy[i][1])

            textcoords = ax.transData.inverted().transform(bbox_expand(curr_bbox, -padding))
            x0 = (textcoords[0][0]+textcoords[1][0]) / 2
            y0 = textcoords[0][1]
            y1 = textcoords[1][1]

            if y1 < pt_data[1]:
                pt_text = (x0, (y1+y0)/2)
            else:
                pt_text = (x0, (y1+y0)/2)

            p = mpatches.FancyArrowPatch(posA=pt_text, posB=pt_data, 
                clip_on=False, **arrowprops)
            ax.add_patch(p)
    
    plt.gcf().canvas.draw()

    # Returns in case more annotations are
    # to be added at a later time
    return bboxes


def get_kwargs(kwargs, keys):
    res = {}
    for k in keys:
        if k in kwargs:
            res[k]=kwargs.pop(k)
    return res

def annotate_span(ax, intervals, loc='bottom', shift=None, padding=5, **kwargs):
    """"""
    x = []
    labels = []

    local_kwargs = get_kwargs(kwargs, ['fc', 'alpha'])

    for interval in intervals['data']:
        mid_pt = (interval.start + interval.end) / 2
        ax.axvspan(interval.start, interval.end, **local_kwargs)

        x.append(mid_pt)
        labels.append(mid_pt)

    if loc not in ['bottom', 'top']:
        raise ValueError("loc must be either 'top' or 'bottom'")
    idx = 0 if loc == 'bottom' else 1

    annotate(ax, labels, list(zip(x, np.repeat(ax.get_ylim()[idx], len(x)))), **kwargs)
