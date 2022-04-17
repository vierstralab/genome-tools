"""General class structure to load and process datasets
"""

from genome_tools.data.loaders import data_loader

import logging
logger = logging.getLogger(__name__)

def iterable_cycle(iterable):
    """
    Args:
      iterable: object with an __iter__ method that can be called multiple times
    """
    while True:
        for x in iterable:
            yield x

class base_dataset(object):
    def batch_iter(self, **kwargs):
        raise  NotImplementedError

    def batch_train_iter(self, cycle=True, **kwargs):
        
        if cycle:
            return ((x["inputs"], x["targets"])
                    for x in iterable_cycle(self._batch_iterable(**kwargs)))
        else:
            return ((x["inputs"], x["targets"]) 
                    for x in self.batch_iter(**kwargs))
        
    def batch_predict_iter(self, cycle=True, **kwargs):
        
        return ((x["inputs"], x["targets"]) 
                    for x in self.batch_iter(**kwargs))
        
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
        
        dp = data_loader(self, 
                        batch_size=batch_size,
                        num_workers=num_workers,
                        **kwargs)

        return dp

    def batch_iter(self, batch_size=1, num_workers=0, **kwargs):
        
        dp =  self._batch_iterable(batch_size=batch_size,
                                    num_workers=num_workers,
                                    **kwargs)
        
        return iter(dp)

    # Example of overloaded `load_all` function
    # def load_all(self, batch_size=1, **kwargs):
    #     """Load all data"""
    #     return numpy_collate_concat([x for x in tqdm(self.batch_iter(batch_size, **kwargs))])
    