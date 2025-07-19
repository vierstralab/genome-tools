import numpy as np
from genome_tools.data.extractors import BigwigExtractor as BigwigExtractor


def smooth_signal(bp_signal, bandwidth, step):
    """
    Smooth signal using a convolution kernel.

    Parameters
    ----------
    bp_signal : np.ndarray
        Per-nucleotide signal.
    bandwidth : int
        Bandwidth for smoothing.
    step : int
        Step for smoothing.

    Returns
    -------
    res : np.ndarray
        Smoothed signal.
    """
    kernel = np.ones(bandwidth) / bandwidth
    smoothed = np.convolve(bp_signal, kernel, mode='same')
    offset = (step - 1) // 2
    res = smoothed[offset::step]
    return res


def smooth_and_aggregate_per_nucleotide_signal(interval, per_nucleotide_signal_files, step, bandwidth):
    """
    Smooth and aggregate per_nucleotide_signal from multiple files.

    Parameters
    ----------
    interval : GenomicInterval
        Genomic interval to extract.
    per_nucleotide_signal_files : list of str
        List of paths to per_nucleotide_signal files.
    step : int
        Step for smoothing.
    bandwidth : int
        Bandwidth for smoothing.
    
    Returns
    -------
    segs : np.ndarray
        Smoothed and averaged signal
    """
    if step is None:
        step = 10
    if bandwidth is None:
        bandwidth = 50
    for i, cc_file in enumerate(per_nucleotide_signal_files):
        segments = BigwigExtractor(cc_file)
        sig = np.nan_to_num(segments[interval])
        if bandwidth > 1:
            sig = smooth_signal(sig, bandwidth, step)
        if i == 0:
            segs = sig
        else:
            segs += sig
    return segs / len(per_nucleotide_signal_files)