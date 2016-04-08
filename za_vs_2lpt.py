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
        while len(args[1]) <3:
            args[1] = '0' + args[1]
        path_2lpt = '/scratch/sissomdj/projects/rockstar/simulations/' + args[0] \
                    + '/2lpt/snap' + args[1] + '/halos/*.bgc2'
        path_za = '/scratch/sissomdj/projects/rockstar/simulations/' + args[0] \
                    + '/za/snap' + args[1] + '/halos/*.bgc2'
        #path_2lpt = '../data/halos_0.1.bgc2'
        #path_2lpt = '../data/halos_0.2.bgc2'
        H2 = HalfMassRadius(path_2lpt, verbose==False)
        HZ = HalfMassRadius(path_za, verbose==False)

        bins = np.linspace(1, 5, 26)    # 25 bins

        H2.read_data()
        H2.filter(100)
        H2.center_halos()
        H2.get_covariance_matrices()
        H2.get_eigenvectors()
        H2.convert_bases()
        H2.get_radii()
        H2.get_half_mass_radii()

        HZ.read_data()
        HZ.filter(100)
        HZ.center_halos()
        HZ.get_covariance_matrices()
        HZ.get_eigenvectors()
        HZ.convert_bases()
        HZ.get_radii()
        HZ.get_half_mass_radii()

        mean2, counts2, std_dev2 = analysis.radius_distribution(H2, 'eval', bins, False )
        meanz, countsz, std_devz = analysis.radius_distribution(HZ, 'eval', bins, False )
        n2 = sum(counts)
        nz = sum(countz)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        _ = ax.hist(mean2, bins, label='2PLT', weights=counts2, color='r', histtype='step')
        _ = ax.hist(meanz, bins, label='ZA', weights=counts2, color='b')
        ax1.set_xlabel('Outer to Inner halo radius')
        ax1.set_ylabel('Normalized Frequency')
        fig.suptitle('2LPT vs. ZA')
        ax1.set_title('2LPT: $\\mu $: ' + '{:.3f}'.format(np.average(mean2, weights=counts2)) + \
                    ', N: ' + str(n2) + '\tZA: $\\mu $: ' + '{:.3f}'.format(np.average(meanz, weights=countsz)) + \
                    ', N: ' + str(nz))
        _ = ax1.legend()
        ax1.grid(True)

        fig.save('results/' + args[0] + '/' + args[1] + '.png')
        return 0


def main(n, box):
    if not os.path.isdir('results'):
        os.makedirs('results')
    if not os.path.isdir('results/' + box):
        os.makedirs('results/' + box)

    M = MyMulticore(n)
    M.args = [(box, str(i)) for i in range(1,2)]
    M.begin()
    _ = M.get_results()
    print 'done'


if __name__=='__main__':
    main(int(sys.argv[1]), sys.argv[2])
