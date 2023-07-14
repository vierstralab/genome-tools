"""Utilities for sequence manipulation"""

import numpy as np

IN_MAP = np.asarray(
    [[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype="int8"
)


def get_base_encoding():
    base_encoding = np.zeros((ord("t") + 1, 1), dtype=np.int)
    base_encoding[ord("A")] = 1
    base_encoding[ord("C")] = 2
    base_encoding[ord("G")] = 3
    base_encoding[ord("T")] = 4
    base_encoding[ord("a")] = 1
    base_encoding[ord("c")] = 2
    base_encoding[ord("g")] = 3
    base_encoding[ord("t")] = 4
    return base_encoding


def one_hot_encode(seqs, revcomp=False):
    base_encoding = get_base_encoding()
    if isinstance(seqs, str) or isinstance(seqs, bytes):
        seqs = [seqs]
    int_seqs = np.concatenate(
        [base_encoding[list(bytearray(str(s), encoding="utf8"))] for s in seqs], axis=1
    )
    out = IN_MAP[int_seqs]
    if revcomp:
        out = np.flip(np.flip(out, axis=0), axis=1)
    return out
