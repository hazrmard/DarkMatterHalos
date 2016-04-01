from ..halos import *
import multiprocessing as mp

class Multicore:
    def __init__(self, processes=1, split_mode='file'):
        self.processes = processes
        self.split_mode = split_mode
        self.pools = []

    def get_cores(self):
        return mp.cpu_count()

    def balance_load(self, files):
        self._setup_pools()
        if self.split_mode=='files':
            H = HalfMassRadius(files)
            H.read_data(level=1)

    def score(self, h, func=lambda x: x.npart**3):
        if type(h) is list:
            score = sum([func(i) for i in h])
        else:
            score = func(h)
        return score

    def _setup_pools(self):
        for _ in range(self.processes):
            self.pools.append(set())
