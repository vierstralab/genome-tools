from genome_tools.plotting.sequence import seq_plot
import numpy as np

def relative_info_content(pwm):
    p = pwm/np.sum(pwm, axis = 1)[:,np.newaxis]
    ic = 2+np.sum(p * np.nan_to_num(np.log2(p)), axis=1)
    ric = p * ic[:, np.newaxis]
    return ric

def plot_motif_logo(motif_pfm, ax=None, offset=-0.5, rc=False, **kwargs):
    pwm = np.loadtxt(motif_pfm)
    if rc:
        pwm = pwm[::-1, ::-1]
    seq_plot(relative_info_content(pwm.T), ax=ax, offset=offset, **kwargs)
    return ax
