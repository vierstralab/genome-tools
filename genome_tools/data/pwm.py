import numpy as np

from genome_tools import VariantInterval

LETTERS = "ACGT"


def read_pfm(file):
    pfm = np.loadtxt(file)
    pfm += 0.001
    pfm /= pfm.sum(axis=0)[np.newaxis,:]
    return pfm


def complement(base):
    _comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return _comp[base]


def reverse_complement(seq):
    _comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # complement each base, then reverse
    return ''.join(_comp[b] for b in reversed(seq))


def relative_info_content(pwm):
    p = pwm / np.sum(pwm, axis = 1)[:, np.newaxis]
    ic = 2 + np.sum(p * np.nan_to_num(np.log2(p)), axis=1)
    ric = p * ic[:, np.newaxis]
    return ric


def seq_logp(mat, seq, bg=None, weights=None):
    """
    Calculate log-probability of a sequence given a PWM matrix.
    Parameters
    ----------
    mat : np.ndarray
        Position weight matrix of shape (4, L).
    seq : str
        Input sequence of length L.
    bg : list or np.ndarray, optional
        Background nucleotide frequencies. If None, uniform background is assumed.
    weights : list or np.ndarray, optional
        Weights for each position in the sequence. If None, all positions are equally weighted.
    
    Returns
    -------
    float
        Log-probability score of the sequence.
    """
    if weights is None:
        weights = np.ones(len(seq))
    if bg is None:
        bg = [0.25, 0.25, 0.25, 0.25]
    res = 0
    for i, c in enumerate(seq):
        j = LETTERS.find(c)
        res += np.log(mat[j,i] / bg[j]) * weights[i]
    return res


def calc_ddg(seq, ref, alt, offset, pfm, **kwargs):
    """
    Calculate delta delta G for a given sequence between ref and alt alleles
    at a given offset using the provided PFM matrix.
    
    Parameters
    ----------
    seq : str
        Input sequence.
    ref : str
        Reference allele.
    alt : str
        Alternate allele.
    offset : int
        Position in the sequence where the alleles differ.
    pfm : np.ndarray
        Position frequency matrix.
    **kwargs
        Additional keyword arguments passed to seq_logp.
    """
    assert len(seq) == pfm.shape[1] 
    assert offset >= 0
    
    ref_seq = list(seq)
    ref_seq[offset] = ref
    ref_score = seq_logp(pfm, ref_seq, **kwargs)
    
    alt_seq = list(seq)
    alt_seq[offset] = alt
    alt_score = seq_logp(pfm, alt_seq, **kwargs)
    return ref_score, alt_score


def get_allelic_scores(
        pfm_matrix: np.ndarray,
        sequence: str,
        ref: str,
        alt: str,
        offset: int,
        orient: str,
    ):
    if orient == '-':
        seq = reverse_complement(sequence)
        ref = complement(ref)
        alt = complement(alt)
        offset = len(sequence) - offset - 1
    else:
        seq = sequence
        ref = ref
        alt = alt
        offset = offset

        # Then return these three columns 
    return calc_ddg(seq, ref, alt, offset, pfm_matrix)