import numpy as np
import itertools

# collections.(Mapping,Sequence) were
# deprecated in newer version of python
try:
    from collections.abc import Mapping, Sequence
except:
    from collections import Mapping, Sequence

# ------------------------

string_classes = (str, bytes)


def _list_collate(stack_fn=list):
    def list_collate_fn(batch):
        if isinstance(batch[0], Mapping):
            return {key: list_collate_fn([d[key] for d in batch]) for key in batch[0]}
        elif isinstance(batch[0], Sequence):
            transposed = zip(*batch)
            return [list_collate_fn(samples) for samples in transposed]
        else:
            return stack_fn(batch)

    return list_collate_fn


list_collate = _list_collate(list)
list_collate.__doc__ = __name__ + ".list_collate"

list_collate_concat = _list_collate(lambda x: list(itertools.chain(*x)))
list_collate_concat.__doc__ = __name__ + ".list_collate_concat"


def _numpy_collate(stack_fn=np.stack):
    def numpy_collate_fn(batch):
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
            return {key: numpy_collate_fn([d[key] for d in batch]) for key in batch[0]}
        elif isinstance(batch[0], Sequence):
            transposed = zip(*batch)
            return [numpy_collate_fn(samples) for samples in transposed]

        raise TypeError(
            (
                "batch must contain tensors, numbers, dicts or lists; found {}".format(
                    type(batch[0])
                )
            )
        )

    return numpy_collate_fn


numpy_collate = _numpy_collate(np.stack)
numpy_collate.__doc__ = __name__ + ".numpy_collate"

numpy_collate_concat = _numpy_collate(np.concatenate)
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
