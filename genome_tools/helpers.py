# Copyright 2019 Jeff Vierstra

import pysam
import gzip
import subprocess
from io import StringIO, TextIOWrapper
import numpy as np
import pandas as pd


magic_dict = {b"\x1f\x8b\x08": "gz", b"\x42\x5a\x68": "bz2", b"\x50\x4b\x03\x04": "zip"}

max_len = max(len(x) for x in magic_dict)


def get_file_type(filename):
    with open(filename, "rb") as f:
        file_start = f.read(max_len)
        for magic, filetype in magic_dict.items():
            if file_start.startswith(magic):
                return filetype
        return None


def open_file(filename):
    file_type = get_file_type(filename)
    if file_type == "gz":
        return gzip.open(filename, mode="rt")
    else:
        return open(filename)


def read_starch(filename, columns=None):
    # Not efficent, try to avoid starch files
    result = subprocess.run(['unstarch', filename], stdout=subprocess.PIPE, text=True, check=True)
    bed_data = pd.read_table(StringIO(result.stdout), header=None)
    ncols = len(bed_data.columns)
    if columns is None:
        columns = ['#chr', 'start', 'end', *np.arange(3, ncols)]
    assert len(columns) == ncols
    bed_data.columns = columns
    return bed_data


def df_to_tabix(df: pd.DataFrame, tabix_path):
    """
    Convert a DataFrame to a tabix-indexed file.
    Renames 'chrom' column to '#chr' if exists.

    Parameters:
        - df: DataFrame to convert to bed format. First columns are expected to be bed-like (chr start end).
        - tabix_path: Path to the tabix-indexed file.

    Returns:
        - None
    """
    with pysam.BGZFile(tabix_path, 'w') as bgzip_out:
        with TextIOWrapper(bgzip_out, encoding='utf-8') as text_out:
            df.rename(columns={'chrom': '#chr'}).to_csv(text_out, sep='\t', index=False)

    pysam.tabix_index(tabix_path, preset='bed', force=True)
