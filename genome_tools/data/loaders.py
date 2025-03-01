"""Enables iterable processing of datasets for multithreading
"""
import sys
import multiprocessing
import threading
import traceback

from genome_tools.data.samplers import (
    random_sampler,
    sequential_sampler,
    minibatch_sampler,
)
from genome_tools.data.utils import list_collate

import logging

logger = logging.getLogger(__name__)

default_collate = list_collate


class exception_wrapper(object):
    """Wraps an exception plus traceback to communicate across threads"""

    def __init__(self, exc_info):
        self.exc_type = exc_info[0]
        self.exc_msg = "".join(traceback.format_exception(*exc_info))


def _worker_loop(process, index_queue, data_queue, collate_fn):
    """Worker loop to load data"""
    while True:
        r = index_queue.get()
        if r is None:
            break
        idx, batch_indices = r

        try:
            batch = collate_fn([process[i] for i in batch_indices])
        except Exception:
            data_queue.put((idx, exception_wrapper(sys.exc_info())))
        else:
            data_queue.put((idx, batch))

    # clean up function that closes filehandlers, etc.
    if hasattr(process, "cleanup"):
        process.cleanup()


class data_loader_iter(object):
    """Iterable data loader class

    Attributes
    ----------


    """

    def __init__(self, loader):
        """


        Parameters
        ----------
        """
        self.process = loader.process
        self.collate_fn = loader.collate_fn
        self.batch_sampler = loader.batch_sampler
        self.num_workers = loader.num_workers
        self.done_event = threading.Event()

        self.sample_iter = iter(self.batch_sampler)

        if self.num_workers > 1:
            logger.info(f"Using {self.num_workers} threads to process data")

            self.index_queue = multiprocessing.SimpleQueue()
            self.data_queue = multiprocessing.SimpleQueue()
            self.batches_outstanding = 0
            self.shutdown = False
            self.send_idx = 0
            self.rcvf_idx = 0
            self.reorder_dict = {}

            self.workers = [
                multiprocessing.Process(
                    target=_worker_loop,
                    args=(
                        self.process,
                        self.index_queue,
                        self.data_queue,
                        self.collate_fn,
                    ),
                )
                for _ in range(self.num_workers)
            ]

            for w in self.workers:
                w.daemon = True
                w.start()

            for _ in range(2 * self.num_workers):
                self._put_indices()
        else:
            logger.info(f"Using a single threads to process data")
            self.shutdown = False

    def __len__(self):
        return len(self.batch_sampler)

    def __next__(self):
        if self.num_workers <= 1:
            indices = next(self.sample_iter)
            batch = self.collate_fn([self.process[i] for i in indices])
            return batch

        if self.rcvf_idx in self.reorder_dict:
            batch = self.reorder_dict.pop(self.rcvf_idx)
            return self._process_next_batch(batch)

        if self.batches_outstanding == 0:
            self._shutdown_workers()
            raise StopIteration

        while True:
            assert not self.shutdown and self.batches_outstanding > 0
            idx, batch = self.data_queue.get()
            self.batches_outstanding -= 1
            if idx != self.rcvf_idx:
                # store out of order samples
                self.reorder_dict[idx] = batch
                continue
            return self._process_next_batch(batch)

    def __iter__(self):
        return self

    def _put_indices(self):
        assert self.batches_outstanding < 2 * self.num_workers
        indices = next(self.sample_iter, None)
        if indices == None:
            return
        self.index_queue.put((self.send_idx, indices))
        self.batches_outstanding += 1
        self.send_idx += 1

    def _process_next_batch(self, batch):
        self.rcvf_idx += 1
        self._put_indices()
        if isinstance(batch, exception_wrapper):
            raise batch.exc_type(batch.exc_msg)
        return batch

    def __getstate__(self):
        raise NotImplementedError("Class cannot be pickled")

    def _shutdown_workers(self):
        if not self.shutdown:
            self.shutdown = True
            self.done_event.set()
            for _ in self.workers:
                self.index_queue.put(None)
                # wait for worker to join?

    def __del__(self):
        if self.num_workers > 1:
            self._shutdown_workers()


class data_loader(object):
    def __init__(
        self,
        process,
        batch_size=1,
        shuffle=False,
        num_workers=0,
        collate_fn=default_collate,
        sampler=None,
        batch_sampler=None,
        drop_last=False,
    ):
        self.process = process
        self.batch_size = batch_size
        self.num_workers = num_workers
        self.collate_fn = collate_fn

        if batch_sampler is not None:
            if batch_size > 1 or shuffle or sampler is not None or drop_last:
                raise ValueError(
                    "batch_sampler is mutually exclusive with batch_size, shuffle, sampler and drop_last"
                )

        if sampler is not None and shuffle:
            raise ValueError("sampler is mutually exclusive with shuffle")

        if batch_sampler is None:
            if sampler is None:
                if shuffle:
                    sampler = random_sampler(process)
                else:
                    sampler = sequential_sampler(process)
            batch_sampler = minibatch_sampler(sampler, batch_size, drop_last)

        self.sampler = sampler
        self.batch_sampler = batch_sampler

        logger.info(
            f"Processing data with batch_size = {self.batch_size:,} resulting in {len(self):,} batches"
        )
        logger.info(f"Using '{self.collate_fn.__doc__}' to collate batch chunks")

    def __iter__(self):
        return data_loader_iter(self)

    def __len__(self):
        return len(self.batch_sampler)

    def on_iter_end(self):
        """
        Callback when one loop through dataset is complete.
        See 'iterable_cycle'. Only used when dataset is in cycle
        mode.
        """
        if hasattr(self.process, "on_iter_end") and callable(
            getattr(self.process, "on_iter_end")
        ):
            self.process.on_iter_end()
        else:
            pass
