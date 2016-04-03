from halos import *
import numpy as np
import matplotlib.pyplot as plt
import time

class MyMulticore(multicore.Multicore):

    def parallel_process(self, halo):
        helpers.do_all(halo)
        return None

    def post_processing(self, H, results):
        bin_edges = np.linspace(0.0, 1.0, 26)
        ratios = [h.half_mass_radius/max(h.radii) for h in H.halos]
        counts, _ = np.histogram(ratios, bin_edges)

        return counts


def main(n):
    H = HalfMassRadius('../data/*.bgc2', verbose=False)
    H.read_data(level=1)
    H.filter(100)
    m = MyMulticore(n)
    m.get_data_from_class(H)
    m.balance_load()
    b = time.clock()
    m.begin()
    res = m.get_results()
    print 'Time passed =', time.clock() - b

    bin_counts = np.array([0] * 25)
    bin_edges = np.linspace(0.0, 1.0, 26)
    bins_mean = [0.5*(bin_edges[i]+bin_edges[i+1]) for i in range(len(bin_edges)-1)]
    for r in res:
        bin_counts += r
    #print 'BIN COUNTS', bin_counts

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    arg_max = np.argmax(bin_counts)
    ax1.set_title('Half Mass Radii - Max:(' + str(bin_counts[arg_max]) + '); Total:' + str(np.sum(bin_counts)))
    _ = ax1.hist(bins_mean, bin_edges, weights=bin_counts)
    ax1.set_xlabel('Inner to Outer halo radius')
    ax1.set_ylabel('Frequency')
    fig.show()
