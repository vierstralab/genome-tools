# Copyright 2019 Jeff Vierstra

import pysam
import gzip
import subprocess
from io import StringIO, TextIOWrapper, BytesIO
import numpy as np
import pandas as pd
import urllib.request
from urllib.parse import urlparse


magic_dict = {b"\x1f\x8b\x08": "gz", b"\x42\x5a\x68": "bz2", b"\x50\x4b\x03\x04": "zip"}

max_len = max(len(x) for x in magic_dict)

_REMOTE_SCHEMES = {"http", "https", "ftp"}


def is_remote(filename):
    return urlparse(str(filename)).scheme in _REMOTE_SCHEMES


def _sniff(head):
    for magic, filetype in magic_dict.items():
        if head.startswith(magic):
            return filetype
    return None

def get_file_type(filename):
    if is_remote(filename):
        with urllib.request.urlopen(filename) as r:
            return _sniff(r.read(max_len))
    with open(filename, "rb") as f:
        return _sniff(f.read(max_len))


def open_file(filename):
    if is_remote(filename):
        # Whole file is buffered in memory: fine for small files (ideograms, chrom sizes),
        # not for large data. Use pysam/tabix for large remote files.
        with urllib.request.urlopen(filename) as r:
            buf = BytesIO(r.read())
        if _sniff(buf.read(max_len)) == "gz":
            buf.seek(0)
            return gzip.open(buf, mode="rt")
        buf.seek(0)
        return TextIOWrapper(buf)

    if get_file_type(filename) == "gz":
        return gzip.open(filename, mode="rt")
    return open(filename)


def read_starch(filename, columns=None):
    # Not efficent, try to avoid starch files. Currently deprecated
    result = subprocess.run(
        ["unstarch", filename], stdout=subprocess.PIPE, text=True, check=True
    )
    bed_data = pd.read_table(StringIO(result.stdout), header=None)
    ncols = len(bed_data.columns)
    if columns is None:
        columns = ["#chr", "start", "end", *np.arange(3, ncols)]
    assert len(columns) == ncols
    bed_data.columns = columns
    return bed_data


def df_to_tabix(df: pd.DataFrame, tabix_path, **kwargs):
    """
    Convert a DataFrame to a tabix-indexed file.
    Renames 'chrom' column to '#chr' if exists.

    Parameters:
        - df: DataFrame to convert to bed format. First columns are expected to be bed-like (chr start end).
        - tabix_path: Path to the tabix-indexed file.
        - **kwargs: Additional keyword arguments to pass to pandas.DataFrame.to_csv

    Returns:
        - None
    """
    with pysam.BGZFile(tabix_path, "w") as bgzip_out:
        with TextIOWrapper(bgzip_out, encoding="utf-8") as text_out:
            df.rename(columns={"chrom": "#chr"}).to_csv(
                text_out, sep="\t", index=False, **kwargs
            )

    pysam.tabix_index(tabix_path, preset="bed", force=True)
