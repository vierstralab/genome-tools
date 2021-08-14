class base_sampler(object):
    def __init__(self, data_source):
        pass

    def __iter__(self):
        raise NotImplementedError

    def __len__(self):
        raise NotImplementedError

class sequential_sampler(base_sampler):
    def __init__(self, data_source):
        self.data_source = data_source
    
    def __iter__(self):
        return iter(range(len(self.data_source)))

    def __len__(self):
        return len(self.data_source)

class minibatch_sampler(base_sampler):
    def __init__(self, sampler, batch_size):
        self.sampler = sampler
        self.batch_size = batch_size

    def __iter__(self):
        batch = []
        for idx in self.sampler:
            batch.append(idx)
            if len(batch) == self.batch_size:
                yield batch
                batch = []
        if len(batch) > 0:
            yield batch
    
    def __len__(self):
        return (len(self.sampler) + self.batch_size - 1) // self.batch_size
        

