from halos import *
import numpy as np
import matplotlib.pyplot as plt
import time
import sys

plt.ioff()

class MyMulticore(multicore.Multicore):

    def parallel_process(self, halo, args=None):
        helpers.do_all(halo)
        return None

    def post_processing(self, H, results, args=None):
        bin_edges = np.linspace(0.0, 1.0, 26)
        ratios = [h.half_mass_radius/max(h.radii) for h in H.halos]
        #ratios = [h.inner_R.x/h.cleave().x for h in H.halos]
        counts, _ = np.histogram(ratios, bin_edges)
        delta_t = time.clock() - self.delta
        return (counts, delta_t)


def main(n):
    b = time.clock()
    H = Halos('../data/Halos_0.1.bgc2', verbose=False)
    #H = Halos('/scratch/sissomdj/projects/simulations/rockstar/box1/za/snap059/halos/*.bgc2', verbose=False)
    H.read_data(level=1)
    H.filter(100)
    m = MyMulticore(n)
    m.score_metric = multicore.quadratic
    m.get_data_from_class(H)
    m.balance_load()
    m.begin()
    res = m.get_results()

    print 'Process start = ', time.strftime('%Y-%m-%d %H:%M:%S')
    print 'Time passed =', '{0:.3f}'.format(time.clock()-b)
    print 'Subprocess #\tHalos\t\tTime /s'
    i=1
    for counts, t in res:
        print str(i) + '\t\t' + str(sum(counts)) + '\t\t' + '{0:.3f}'.format(t)
        i+=1
    print
    print 'Source files:'
    for f in m.files:
        print f
    print

    bin_counts = np.array([0] * 25)
    bin_edges = np.linspace(0.0, 1.0, 26)
    bins_mean = [0.5*(bin_edges[i]+bin_edges[i+1]) for i in range(len(bin_edges)-1)]
    for r, _ in res:
        bin_counts += r

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    arg_max = np.argmax(bin_counts)
    ax1.set_title('Half Mass Radii - Max:(' + str(bin_counts[arg_max]) + '); Total:' + str(np.sum(bin_counts)))
    _ = ax1.hist(bins_mean, bin_edges, weights=bin_counts)
    ax1.set_xlabel('Inner to Outer halo radius')
    ax1.set_ylabel('Frequency')
    fig.savefig('result.png')


if __name__=='__main__':
    #num_processes = MyMulticore().get_cores()
    num_processes=int(sys.argv[1])
    print 'Number of processes:', num_processes
    main(num_processes)
    print 'Finished'
