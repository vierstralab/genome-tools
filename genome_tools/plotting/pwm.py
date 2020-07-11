from matplotlib.patheffects import RendererBase
from matplotlib.transforms import offset_copy
from matplotlib.patches import Rectangle

import numpy as np

import seaborn as sns

class Scale(RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy
        
    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy) + affine
        renderer.draw_path(gc, tpath, affine, rgbFace)
        
basic = "ACGT"
iupac = "XACMGRSVTWYHKDBN"

bg_uni =  np.array([0.25, 0.25, 0.25, 0.25])
bg = np.array([0.2, 0.3, 0.3, 0.2])

iupac_colors = ['black'] * 16
iupac_colors[1<<0] = 'green'
iupac_colors[1<<1] = 'blue'
iupac_colors[1<<2] = 'gold'
iupac_colors[1<<3] = 'red'


class pwm(object):

    def __init__(self, mat):
        self.data = mat

    def rel_info_content(self, bg=None):
        mat = np.copy(self.data).T
        p=mat/np.sum(mat, axis = 1)[:,np.newaxis]
        if bg is None:
            ic=2+np.sum(p*np.nan_to_num(np.log2(p)), axis = 1)
        else:
            ic=2+np.sum(p*np.nan_to_num(np.log2(p/bg[np.newaxis,:])), axis = 1)
            
        ric=p*ic[:,np.newaxis]

        return ric

    def seq_logp(self, seq, bg):
        res=0
        for i, c in enumerate(seq):
            j = basic.find(c)
            res += np.log(self.data[i, j]/bg[j])
        return res

    def render(self, fig, ax, pad=0, xoffset=0, xlim=None, type='default', bg=None, rc=False):

        if type=='ic':
            mat=self.rel_info_content(bg=bg).T
            ax.set_ylim(0, 2.1)
        elif type=='affinity':
            ddg = np.log(self.data+1e-3)
            mat=ddg-np.mean(ddg, axis=0)[np.newaxis,:]
            top=np.max(np.ma.sum(np.ma.array(mat,  mask=mat<0), axis=0))
            bottom=np.min(np.ma.sum(np.ma.array(mat,  mask=mat>=0), axis=0))
            ax.set_ylim(bottom, top)
        else:
            mat=self.data
            ax.set_ylim(0, 1.05)
        
        if rc:
            mat = mat[::-1,:][:,::-1]
            

        w = mat.shape[1]
        if xlim:
            ax.set_xlim(xlim[0], xlim[1])
        else:
            ax.set_xlim(xoffset-pad, xoffset+w+pad)
        
        wscale=1
        init = False

        for i in range(w):

            heights=mat[:,i]
            order=np.argsort(heights)
            reorder=np.concatenate([order[heights[order]>0], order[heights[order]<=0]])

            yshift_pos = 0
            yshift_neg = np.sum(heights[heights<=0])

            for j, b in enumerate(reorder):
                base=iupac[(1<<b)]
                color=iupac_colors[(1<<b)]
                scale=heights[b]

                if scale>0:
                    yshift=yshift_pos
                    alpha=1
                else:
                    yshift=yshift_neg
                    alpha=0.4

                h=np.abs(scale)

                t = ax.text(i+xoffset, yshift, base, ha='left', va='baseline', 
                    color=color, fontsize=40, family='monospace', 
                    fontname="Bitstream Vera Sans Mono", weight='bold', alpha=alpha)

                if not init:
                    fig.canvas.draw()
                    ext = t.get_window_extent(t._renderer)

                    x0, x1, y0, y1 = ext.x0, ext.x1, ext.y0, ext.y1
                    nx0, ny0 = ax.transData.inverted().transform((x0, y0))
                    nx1, ny1 = ax.transData.inverted().transform((x1, y1))
                    
                    wscale=1.0/(nx1-nx0)
                    hscale=(1.0/(ny1))
                    
                    height=ext.height        
                    
                    init=True

                t.set_path_effects([Scale(wscale, hscale*h)])
                
                if scale>0:
                    yshift_pos+=h
                else:
                    yshift_neg+=h
                
        ax.set_xticks(np.arange(xoffset, xoffset+mat.shape[1])+0.5)
        ax.set_xticklabels([])
        ax.xaxis.set_visible(False)
        
        ax.set_yticks([])
        ax.set_yticklabels([])
        #ax.yaxis.set_visible(False)
        
        sns.despine(ax=ax, left=True, bottom=True)

    