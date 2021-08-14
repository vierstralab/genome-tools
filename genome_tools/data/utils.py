import numpy as np
import collections
import itertools

string_classes = (str, bytes)

def _list_collate(stack_fn=list):
    def list_collate_fn(batch):
        if isinstance(batch[0], collections.Mapping):
            return {key: list_collate_fn([d[key] for d in batch]) for key in batch[0]}
        elif isinstance(batch[0], collections.Sequence):
            transposed = zip(*batch)
            return [list_collate_fn(samples) for samples in transposed]
        else:
            return stack_fn(batch)

    return list_collate_fn

list_collate = _list_collate(list)
list_collate.__doc__ = __name__ + '.list_collate'

list_collate_concat = _list_collate(lambda x: list(itertools.chain(*x)))
list_collate_concat.__doc__ = __name__ + '.list_collate_concat'

def _numpy_collate(stack_fn=np.stack):
    def numpy_collate_fn(batch):
        "Puts each data field into a tensor with outer dimension batch size"
        if type(batch[0]).__module__ == 'numpy':
            elem = batch[0]
            if type(elem).__name__ == 'ndarray':
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
        elif isinstance(batch[0], collections.Mapping):
            return {key: numpy_collate_fn([d[key] for d in batch]) for key in batch[0]}
        elif isinstance(batch[0], collections.Sequence):
            transposed = zip(*batch)
            return [numpy_collate_fn(samples) for samples in transposed]

        raise TypeError(("batch must contain tensors, numbers, dicts or lists; found {}"
                         .format(type(batch[0]))))

    return numpy_collate_fn

numpy_collate = _numpy_collate(np.stack)
numpy_collate.__doc__ = __name__ + '.numpy_collate'

numpy_collate_concat = _numpy_collate(np.concatenate)
numpy_collate_concat.__doc__ = __name__ + '.numpy_collate_concat'