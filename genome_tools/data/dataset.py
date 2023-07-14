"""General class structure to load and process datasets
"""

from genome_tools.data.loaders import data_loader
from genome_tools.data.utils import (
    iterable_cycle,
    numpy_collate,
    numpy_collate_concat,
    get_dataset_lens,
    get_dataset_item,
)

from tqdm import tqdm

import logging

logger = logging.getLogger(__name__)


class base_dataset(object):
    def batch_iter(self, **kwargs):
        raise NotImplementedError

    def batch_train_iter(self, cycle=True, **kwargs):
        if cycle:
            return (
                (x["inputs"], x["targets"], x["weights"])
                for x in iterable_cycle(self._batch_iterable(**kwargs))
            )
        else:
            return (
                (x["inputs"], x["targets"], x["weights"])
                for x in self.batch_iter(**kwargs)
            )

    def batch_predict_iter(self, **kwargs):
        return (x["inputs"] for x in self.batch_iter(**kwargs))

    def batch_eval_iter(self, **kwargs):
        return (
            (x["inputs"], x["targets"], x["metadata"])
            for x in self.batch_iter(**kwargs)
        )

    def load_all(self, **kwargs):
        """Load an entire dataset"""
        return [x for x in self.batch_iter(**kwargs)]


class dataset(base_dataset):
    """All datasets should subclass this class. All subclases should
    override `__len__` and `__getitem__` to support integer indexing"""

    def __getitem__(self, index):
        """Return one sample"""
        raise NotImplementedError

    def __len__(self):
        """Return number of all elements"""
        raise NotImplementedError

    def _batch_iterable(self, batch_size=1, num_workers=0, **kwargs):
        dl = data_loader(self, batch_size=batch_size, num_workers=num_workers, **kwargs)

        return dl

    def batch_iter(self, batch_size=1, num_workers=0, **kwargs):
        dl = self._batch_iterable(
            batch_size=batch_size, num_workers=num_workers, **kwargs
        )

        return iter(dl)

    # Example of overloaded `load_all` function
    def load_all(self, batch_size=1, **kwargs):
        """Load all data"""
        return numpy_collate_concat(
            [x for x in tqdm(self.batch_iter(batch_size, **kwargs))]
        )


class preloaded_dataset(base_dataset):
    data_fn = None

    @classmethod
    def from_fn(cls, data_fn):
        cls.data_fn = staticmethod(data_fn)
        return cls

    @classmethod
    def from_data(cls, data):
        return cls.from_data_fn(lambda: data)()

    @classmethod
    def _get_data_fn(cls):
        assert cls.data_fn is not None
        return cls.data_fn

    def __init__(self, *args, **kwargs):
        self.data = self._get_data_fn()(*args, **kwargs)
        lens = get_dataset_lens(self.data, require_numpy=True)
        assert len(set(lens)) == 1
        self.n = lens[0]

    def __len__(self):
        return self.n

    def __getitem__(self, index):
        return get_dataset_item(self.data, index)

    def _batch_iterable(self, batch_size=32, shuffle=False, drop_last=False, **kwargs):
        dl = data_loader(
            self,
            batch_size=batch_size,
            collate_fn=numpy_collate,
            shuffle=shuffle,
            num_workers=0,
            drop_last=drop_last,
        )
        return dl

    def batch_iter(self, batch_size=32, shuffle=False, drop_last=False, **kwargs):
        dl = self._batch_iterable(
            batch_size=batch_size, shuffle=shuffle, drop_last=drop_last, **kwargs
        )
        return iter(dl)

    def load_all(self, **kwargs):
        return self.data
