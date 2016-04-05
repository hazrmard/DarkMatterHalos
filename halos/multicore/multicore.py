from ..halos import *
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import os
import time
plt.ioff()

class Multicore:
    def __init__(self, processes=1):
        self.processes = processes      # number of processes to spawn
        self.pools = []                 # list of lists of halos/process
        self.h = []                     # list of halos prior to balancing
        self.score_metric = cubic       # determines score of each halo
        self.files = None                   # source files
        self.queue = mp.Queue()             # inter process communication
        self.process_list = []              # list of spawned processes
        self.delta = 0                      # time stamp of worker start

    def get_cores(self):
        """Returns number of cores (including virtual) on the machine.
        """
        return mp.cpu_count()

    def get_data_from_files(self, files):
        """Loads halo metadata from files directly.
        :files single or list of paths. Paths can have wildcards.
        """
        H = HalfMassRadius(files, verbose=False)
        H.read_data(level=1)
        self.h = H.h
        self.files = H.files

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

    def visualize(self, fpath=None):
        """Show a bar graph containing halos/process and load/process
        :fpath output path of saved figure
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
        if fpath is None:
            fig.show()
        else:
            fig.savefig(fpath)

    def worker(self, pool_ids, files, queue):
        """The target function of each process. Runs parallel_process() over
        each halo in thread and passes the results to post_processing() for
        aggregation etc. Final result as determined by post_processing() is given
        to parent process.
        :pool_ids a set of halo ids to be loaded to current process
        :files file paths to all source data
        :queue a multiprocessing.Queue() instance. For inter-process communication
        """
        self.delta = time.clock()
        H = HalfMassRadius(files, verbose=False)
        H.read_data(level=2, sieve=pool_ids)
        results = []
        for halo in H.halos:
            results.append(self.parallel_process(halo))
        final_result = self.post_processing(H, results)
        queue.put(final_result)

    def parallel_process(self, halo):
        """the function applied to a single Halo instance. This function is called
        in a loop for all halos in a thread. The return value is appended to a list of
        results of parallel_process() on all halos in a thread.
        :halo the Halo instance automaticall passed to this function
        """
        return None

    def post_processing(self, H, results):
        """The HalfMassRadius object of a particular thread and the list of results
        of parallel_process() on all halos in a thread are passed to this function.
        can be used for aggregating values or making final edits to results before being
        passed on to the main process. The return value is appended to a list of results
        for each process on the main process.
        """
        return None

    def begin(self):
        """spawns processes. Each processes instantiates a HalfMassRadius object
        containing halos as determined by the balance_load() function.
        """
        for i in range(self.processes):
            pool_ids = set([halo.id for halo in self.pools[i]])
            if len(pool_ids)==0:
                pool_ids = None
            p = mp.Process(target=self.__class__(1).worker, args=(pool_ids, self.files, self.queue))
            self.process_list.append(p)
            p.start()

    def get_results(self):
        """wait for processes to finish and obtain the results of post_processing()
        on each process in a list.
        """
        results = []
        for _ in range(self.processes):
            results.append(self.queue.get())
        for p in self.process_list:
            p.terminate()
            p.join()
        return results

    def score(self, h, score_metric):
        """a function used to compute the cost of a particular halo/collection
        of halos. This is used to balance loads across each process.
        :h a Halo object on a collection of Halo objects
        :score_metric a function that accepts a Halo object and returns its cost.
        Can be changed permanently by editing Multicore.score_metric attribute.
        """
        #func = self.score_metric if self.score_metric is not None else func
        if hasattr(h, '__iter__'):
            score = sum([score_metric(i) for i in h])
        else:
            score = score_metric(h)
        return score

    def _setup_pools(self):
        """creates a list of lists. Each sublist contains metadata of halos going
        in each process.
        """
        self.pools = []
        for _ in range(self.processes):
            self.pools.append([])

    def _sort_halos(self, score_metric):
        """sorts halo objects with metadata in ascending order of cost.
        """
        sorted_indices = np.argsort([self.score(halo, score_metric) for halo in self.h])
        self.h = list(np.array(self.h)[sorted_indices])


def cubic(x):
    return x.npart**3


def quadratic(x):
    return x.npart**2


def linear(x):
    return x.npart
