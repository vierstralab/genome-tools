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


def seq_logp(mat, seq, bg=None):
    if bg is None:
        bg = [0.25, 0.25, 0.25, 0.25]
    res = 0
    for i, c in enumerate(seq):
        j = LETTERS.find(c)
        res += np.log(mat[j,i]/bg[j])
    return res


def calc_ddg(seq, ref, alt, offset, pfm):
    assert len(seq) == pfm.shape[1] 
    assert offset >= 0
    
    ref_seq = list(seq)
    ref_seq[offset] = ref
    ref_score = seq_logp(pfm, ref_seq)
    
    alt_seq = list(seq)
    alt_seq[offset] = alt
    alt_score = seq_logp(pfm, alt_seq)
    return ref_score, alt_score

def get_allelic_scores(
        pfm_matrix: np.ndarray,
        motif_start,
        motif_end,
        sequence: str,
        orient,
        variant_interval: VariantInterval
    ):
    # constants from cavs
    assert motif_end - motif_start == pfm_matrix.shape[1]
    if orient == '-':
        seq = reverse_complement(sequence)
        ref = complement(variant_interval.ref)
        alt = complement(variant_interval.alt)
        offset = motif_end - variant_interval.end
    else:
        seq = sequence
        ref = variant_interval.ref
        alt = variant_interval.alt
        offset = variant_interval.start - motif_start

        # Then return these three columns 
    return calc_ddg(seq, ref, alt, offset, pfm_matrix)