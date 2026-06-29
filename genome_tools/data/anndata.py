import zarr
import anndata as ad
import dask.array as da
import numpy as np
import scipy.sparse as sp


def _read_zarr_group(dat):
    if isinstance(dat, zarr.Array):
        return da.from_zarr(dat)
    else:
        return ad.experimental.sparse_dataset(dat).to_memory()


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
