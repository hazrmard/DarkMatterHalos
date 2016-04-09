from halos import HalfMassRadius
from halos.multicore import Multicore
import analysis
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.ioff()

class MyMulticore(Multicore):

    def worker(self, pool_ids, files ,queue, args):
        try:
            while len(args[1]) <3:
                args[1] = '0' + args[1]
            path_2lpt = '/scratch/sissomdj/projects/simulations/rockstar/' + args[0] \
                        + '/2lpt/snap' + args[1] + '/halos/*1.bgc2'
            path_za = '/scratch/sissomdj/projects/simulations/rockstar/' + args[0] \
                        + '/za/snap' + args[1] + '/halos/*1.bgc2'
            #path_2lpt = '../data/halos_0.1.bgc2'
            #path_2lpt = '../data/halos_0.2.bgc2'
            H2 = HalfMassRadius(path_2lpt, verbose=False)
            HZ = HalfMassRadius(path_za, verbose=False)

            bins = np.linspace(1, 5, 26)    # 25 bins
            bins_mean = [0.5*(bins[i]+bins[i+1]) for i in range(len(bins)-1)]

            H2.read_data()
            H2.filter(100)
            H2.center_halos()
            H2.get_covariance_matrices()
            H2.get_eigenvectors()
            H2.convert_bases()
            H2.get_radii()
            H2.get_half_mass_radii()
            ratios2 = [max(h.radii)/h.half_mass_radius for h in H2.halos]
            del H2

            HZ.read_data()
            HZ.filter(100)
            HZ.center_halos()
            HZ.get_covariance_matrices()
            HZ.get_eigenvectors()
            HZ.convert_bases()
            HZ.get_radii()
            HZ.get_half_mass_radii()
            ratiosz = [max(h.radii)/h.half_mass_radius for h in HZ.halos]
            del HZ

            counts2, _ = np.histogram(ratios2, bins)
            countsz, _ = np.histogram(ratiosz, bins)
            n2 = sum(counts2)
            nz = sum(countsz)

            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            _ = ax1.hist(bins_mean, bins, label='2LPT', weights=counts2, color='r', histtype='step')
            _ = ax1.hist(bins_mean, bins, label='ZA', color='b',weights=countsz, histtype='step')
            ax1.set_xlabel('Outer to Inner halo radius')
            ax1.set_ylabel('Frequency')
            fig.suptitle('2LPT vs. ZA')
            ax1.set_title('2LPT: $\\mu $= ' + '{:.3f}'.format(np.mean(ratios2)) + \
                        ', N= ' + str(n2) + '\tZA: $\\mu $= ' + '{:.3f}'.format(np.mean(ratiosz)) + \
                        ', N= ' + str(nz))
            _ = ax1.legend()
            ax1.grid(True)

            fig.savefig('results/' + args[0] + '/' + args[1] + '.png')
            queue.put((0,None))
            return (0, None)
        except Exception as e:
            queue.put((-1, e))
            return -1


def main(n, box):
    if not os.path.isdir('results'):
        os.makedirs('results')
    if not os.path.isdir('results/' + box):
        os.makedirs('results/' + box)

    print '\nRunning', n, 'processes on', box, '\n'
    M = MyMulticore(n)
    M.args = [[box, str(i)] for i in range(0,62)]
    M.begin()
    status = M.get_results()
    for s, msg in status:
        if s==-1:
            print e
            print '\n\n'
    print 'done'


if __name__=='__main__':
    main(int(sys.argv[1]), sys.argv[2])
