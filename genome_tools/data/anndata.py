import zarr
import anndata as ad
import dask.array as da
import numpy as np
import scipy.sparse as sp


def _read_zarr_group(dat):
    if type(dat) == zarr.Array:
        return da.from_zarr(dat)
    else:
        return ad.experimental.sparse_dataset(dat).to_memory()


def read_zarr_backed(filename):
    d = {}
    attributes = ["obsm", "varm", "obsp", "varp", "uns"]
    df_attributes = ["obs", "var"]

    with zarr.open(filename, "r") as f:
        if "encoding-type" in f.attrs:
            attributes.extend(df_attributes)

        if 'X' in f:
            d['X'] = _read_zarr_group(f['X'])

        for k in attributes:
            if k not in f:
                continue
            d[k] = ad.experimental.read_elem(f[k])

        if 'layers' in f.keys():
            d['layers'] = {}
            for layer in f['layers']:
                d['layers'][layer] = _read_zarr_group(f['layers'][layer])

        adata = ad.AnnData(**d)

    return adata
