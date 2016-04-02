from ..halos import *
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import os

class Multicore:
    def __init__(self, processes=1):
        self.processes = processes
        self.pools = []
        self.h = []
        self.score_metric = default_score
        self.files = None
        self.queue = mp.Queue()
        self.process_list = []

    def get_cores(self):
        """Returns number of cores (including virtual) on the machine.
        """
        return mp.cpu_count()

    def get_data_from_files(self, files):
        """Loads halo metadata from files directly.
        :files single or list of paths. Paths can have wildcards.
        """
        H = HalfMassRadius(files)
        H.read_data(level=1)
        self.h = H.h
        self.files = files

    def get_data_from_class(self, c):
        """Extract halo metadeta from a HalfMassRadius object.
        :c a HalfMassRadius instance
        """
        self.h = c.h
        self.files = c.files

    def balance_load(self, score_metric=None):
        """given the number of processes to run, assigns halos to each process
        so as to even out load per process. Load is determined by a score metric.
        Default score metric is O(n^3) where n is # of particles in halo.
        """
        score_metric = self.score_metric if score_metric is None else score_metric
        self._setup_pools()
        self._sort_halos(score_metric)
        pool_scores = [self.score(pool, score_metric) for pool in self.pools]
        min_pool = self.pools[0]
        min_pool_index = 0
        while self.h:
            halo = self.h.pop()
            min_pool.append(halo)
            pool_scores[min_pool_index] += self.score(halo, score_metric)
            min_pool_index = np.argmin(pool_scores)
            min_pool = self.pools[min_pool_index]
        #print 'halos/pool', [len(pool, score_metric) for pool in self.pools]
        #print 'score/pool', [self.score(pool, score_metric) for pool in self.pools]

    def visualize(self):
        """Show a bar graph containing halos/process and load/process
        """
        d1 = [len(pool) for pool in self.pools]
        d2 = [self.score(pool, self.score_metric) for pool in self.pools]
        bar_width = 0.35
        index = np.arange(len(d1))
        fig, ax1 = plt.subplots()
        ax1.bar(index, d1, bar_width, label='Halos in Pool', color='b')
        ax1.set_ylabel('Halos in Pool', color='b')
        for tl in ax1.get_yticklabels():
            tl.set_color('b')
        ax1.set_xticks(index+bar_width)
        ax1.set_xticklabels([str(x) for x in range(len(self.pools))])
        ax1.set_title('Load Distribution in Pools')
        ax1.set_xlabel('Pool #')

        ax2 = ax1.twinx()
        ax2.bar(index+bar_width, d2, bar_width, label='Pool Load', color='r')
        ax2.set_ylabel('Pool Load', color='r')
        for tl in ax2.get_yticklabels():
            tl.set_color('r')

        fig.show()

    def worker(self, pool_ids, files, queue):
        H = HalfMassRadius(files)
        #lock.acquire()
        #print os.getpid(), "process reading data"
        H.read_data(level=2, sieve=pool_ids)
        #print os.getpid()
        #lock.release()
        results = []
        for halo in H.halos:
            results.append(self.parallel_process(halo))
        final_result = self.post_processing(H, results)
        queue.put(final_result)

    def parallel_process(self, halo):
        return 1

    def post_processing(self, H, results):
        return sum(results)

    def begin(self):
        for i in range(self.processes):
            pool_ids = set([halo.id for halo in self.pools[i]])
            if len(pool_ids)==0:
                pool_ids = None
            p = mp.Process(target=self.__class__(1).worker, args=(pool_ids, self.files, self.queue))
            self.process_list.append(p)
            p.start()

    def get_results(self):
        results = []
        for _ in range(self.processes):
            results.append(self.queue.get())
        for p in self.process_list:
            p.join()
        return results

    def score(self, h, score_metric):
        #func = self.score_metric if self.score_metric is not None else func
        if hasattr(h, '__iter__'):
            score = sum([score_metric(i) for i in h])
        else:
            score = score_metric(h)
        return score

    def _setup_pools(self):
        self.pools = []
        for _ in range(self.processes):
            self.pools.append([])

    def _sort_halos(self, score_metric):
        sorted_indices = np.argsort([self.score(halo, score_metric) for halo in self.h])
        self.h = list(np.array(self.h)[sorted_indices])


def default_score(x):
    return x.npart**3
