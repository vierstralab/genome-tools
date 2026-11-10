import numpy as np
from genome_tools import GenomicInterval

def remap_columns(matrix: np.ndarray, src_cols: np.ndarray) -> np.ndarray:
    """Output position i = source position src_cols[i] along the last axis; gaps stay zero."""
    valid = (src_cols >= 0) & (src_cols < matrix.shape[-1])
    aligned = np.zeros(matrix.shape[:-1] + (len(src_cols),), dtype=matrix.dtype)
    aligned[..., valid] = matrix[..., src_cols[valid]]
    return aligned


def realign_matrix(
        matrix: np.ndarray,
        matrix_interval: GenomicInterval,
        target_interval: GenomicInterval,
    ) -> np.ndarray:
    """Same-species: source column is a fixed offset from target position."""
    assert matrix.shape[-1] == len(matrix_interval)
    offset = target_interval.start - matrix_interval.start
    src_cols = np.arange(len(target_interval)) + offset

    n_valid = int(((src_cols >= 0) & (src_cols < matrix.shape[-1])).sum())
    n_target = len(target_interval)
    if n_valid == 0:
        print(f"Warning: {target_interval} does not overlap {matrix_interval}; returning zeros.")
    elif n_valid < n_target:
        print(f"Warning: partial overlap ({n_valid}/{n_target}) between {target_interval} and {matrix_interval}; padding with zeros.")

    return remap_columns(matrix, src_cols)



### DEFUNC
# ------------------------
try:
    from collections.abc import Mapping, Sequence
except ImportError:
    from collections import Mapping, Sequence

import functools
import itertools
# collections.(Mapping,Sequence) were
# deprecated in newer version of python

string_classes = (str, bytes)


# Collate functions are module-level functions bound with functools.partial (not closures),
# so they can be pickled for multiprocessing with the spawn start method (macOS default)
def _list_collate_fn(batch, stack_fn=list):
    if isinstance(batch[0], Mapping):
        return {key: _list_collate_fn([d[key] for d in batch], stack_fn) for key in batch[0]}
    elif isinstance(batch[0], Sequence):
        transposed = zip(*batch)
        return [_list_collate_fn(samples, stack_fn) for samples in transposed]
    else:
        return stack_fn(batch)


def _chain_list(x):
    return list(itertools.chain(*x))


list_collate = functools.partial(_list_collate_fn, stack_fn=list)
list_collate.__doc__ = __name__ + ".list_collate"

list_collate_concat = functools.partial(_list_collate_fn, stack_fn=_chain_list)
list_collate_concat.__doc__ = __name__ + ".list_collate_concat"


def _numpy_collate_fn(batch, stack_fn=np.stack):
        "Puts each data field into a tensor with outer dimension batch size"
        if type(batch[0]).__module__ == "numpy":
            elem = batch[0]
            if type(elem).__name__ == "ndarray":
                return stack_fn(batch, 0)
            if elem.shape == ():  # scalars
                return np.array(batch)
        elif isinstance(batch[0], int):
            return np.asarray(batch)
        elif isinstance(batch[0], float):
            return np.asarray(batch)
        elif batch[0] is None:
            return np.asarray(batch)
        elif isinstance(batch[0], string_classes):
            # Also convert to a numpy array
            return np.asarray(batch)
            # return batch
        elif isinstance(batch[0], Mapping):
            return {key: _numpy_collate_fn([d[key] for d in batch], stack_fn) for key in batch[0]}
        elif isinstance(batch[0], Sequence):
            transposed = zip(*batch)
            return [_numpy_collate_fn(samples, stack_fn) for samples in transposed]

        raise TypeError(
            (
                "batch must contain tensors, numbers, dicts or lists; found {}".format(
                    type(batch[0])
                )
            )
        )


numpy_collate = functools.partial(_numpy_collate_fn, stack_fn=np.stack)
numpy_collate.__doc__ = __name__ + ".numpy_collate"

numpy_collate_concat = functools.partial(_numpy_collate_fn, stack_fn=np.concatenate)
numpy_collate_concat.__doc__ = __name__ + ".numpy_collate_concat"

# ------------------------


def check_input_dimensions(intervals):
    """Checks that each interval is of the same width"""
    sizes = list(map(len, intervals))
    if len(np.unique(sizes)) != 1:
        raise IndexError("Inputs are not of same length!")


def iterable_cycle(iterable):
    """
    Args:
      iterable: object with an __iter__ method that can be called multiple times
    """
    while True:
        for x in iterable:
            yield x

        # Run the callback (if exists) after each pass
        # through the dataset
        if hasattr(x, "on_iter_end") and callable(getattr(x, "on_iter_end")):
            x.on_iter_end()


def get_dataset_lens(data, require_numpy=False):
    if type(data).__module__ == "numpy":
        if require_numpy and not data.shape:
            raise ValueError("all numpy arrays need to have at least one axis")
        return [1] if not data.shape else [data.shape[0]]
    elif isinstance(data, int) and not require_numpy:
        return [1]
    elif isinstance(data, float) and not require_numpy:
        return [1]
    elif isinstance(data, string_classes) and not require_numpy:
        # Also convert to a numpy array
        return [1]
        # return data
    elif isinstance(data, Mapping) and not type(data).__module__ == "numpy":
        return sum([get_dataset_lens(data[key], require_numpy) for key in data], [])
    elif isinstance(data, Sequence) and not type(data).__module__ == "numpy":
        return sum([get_dataset_lens(sample, require_numpy) for sample in data], [])
    else:
        raise ValueError("Leafs of the nested structure need to be numpy arrays")


def get_dataset_item(data, idx):
    if type(data).__module__ == "numpy":
        return data[idx]
    elif isinstance(data, Mapping):
        return {key: get_dataset_item(data[key], idx) for key in data}
    elif isinstance(data, Sequence):
        return [get_dataset_item(sample, idx) for sample in data]
    else:
        raise ValueError("Leafs of the nested structure need to be numpy arrays")
