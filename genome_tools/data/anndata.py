from typing import List

import zarr
import anndata as ad
import dask.array as da
import numpy as np
import scipy.sparse as sp


def _read_zarr_group(dat):

    if isinstance(dat, zarr.Array):
        return da.from_zarr(dat)
    elif isinstance(dat, zarr.Group):
        enc = dat.attrs.get("encoding-type", "")
        if enc == "dataframe":
            return ad.experimental.read_elem(dat)  # returns pd.DataFrame
        else:
            return ad.experimental.sparse_dataset(dat).to_memory()
    else:
        raise ValueError('Not a zarr array or zarr group')


def read_zarr_backed(filename, lazy_attrs=('layers', 'varm')):
    d = {}
    attributes = ["layers", "varm", "obsm", "obsp", "varp", "uns"]
    df_attributes = ["obs", "var"]

    with zarr.open(filename, "r") as f:
        if "encoding-type" in f.attrs:
            attributes.extend(df_attributes)

        if 'X' in f:
            d['X'] = _read_zarr_group(f['X'])

        for ad_attr in attributes:
            if ad_attr not in f:
                continue
            elif ad_attr in lazy_attrs:
                d[ad_attr] = {}
                for layer in f[ad_attr]:
                    d[ad_attr][layer] = _read_zarr_group(f[ad_attr][layer])
            else:
                d[ad_attr] = ad.experimental.read_elem(f[ad_attr])

        adata = ad.AnnData(**d)

    return adata


def convert_to_sparse_if_sufficiently_sparse(matrix, threshold=0.7):
    """
    Convert a matrix to sparse format if the proportion of zero values exceeds the threshold.
    """
    if sp.issparse(matrix):
        return matrix

    sparsity = np.count_nonzero(matrix == 0) / matrix.size
    return sp.csr_matrix(matrix) if sparsity > threshold else matrix


def get_matrices_mapping_by_types(matrices: List[str], matrices_types):
    """
    Map each matrix file path to the matrix_type (from matrices_types) its filename encodes,
    based on the '<matrix_type>.npy' / '<matrix_type>.matrix.npy' naming convention.
    """
    matrices_mapping = {}
    for matrix in matrices:
        for matrix_type in matrices_types:
            if matrix.endswith(f'.{matrix_type}.npy') or matrix.endswith(f'{matrix_type}.matrix.npy'):
                matrices_mapping[matrix_type] = matrix
                break
        else:
            raise ValueError(f"Matrix {matrix} not recognized. Expected one of {matrices_types}")
    return matrices_mapping


def add_matrices_to_anndata(adata, matrices_mapping: dict, mask=None):
    """
    Load each matrix file in matrices_mapping and add it as a layer to adata, keyed by its mapping key.
    Each matrix file is expected to be var x obs on disk (transposed to obs x var to match adata).
    """
    for matrix_name, matrix in matrices_mapping.items():
        if matrix_name in adata.layers:
            print(f"Layer {matrix_name} already exists in anndata. Skipping.")
            continue
        matrix = np.load(matrix).T
        if matrix.dtype == np.float64:
            matrix = matrix.astype(np.float32)

        if matrix.shape[0] != adata.shape[0]:
            raise ValueError(f"Matrix {matrix_name} has {matrix.shape[0]} samples, but the number of samples in anndata is {adata.shape[0]}")
        if mask is not None:
            assert mask.sum() == matrix.shape[1]
            data = np.empty(adata.shape, dtype=np.float32)
            data[:] = np.nan
            data[:, mask] = matrix
            matrix = data

        adata.layers[matrix_name] = convert_to_sparse_if_sufficiently_sparse(matrix)


def get_mask_from_column_name(adata, mask_name):
    if mask_name != "":
        mask = adata.var[mask_name]
    else:
        mask = np.ones(adata.var.shape[0], dtype=bool)
    return mask
