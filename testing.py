__author__ = 'Ibrahim'

# This file contains examples of functions that can be used to process halo data.

from halos import *
from halos.helpers import *
import halos.multicore as m
import time
import multiprocessing as mp
import os
import random


def bgc2_test(path='..\data\halos_0.1.bgc2'):
    """
    reading data from a sample bgc2 file containing multiple halos
    :param path: path to file
    :return: a HalfMassRadius instance containing list of halos (<return_variable>.halos)
    """
    t = HalfMassRadius(path)
    t.read_data()
    t.filter(100)         # filter out halos w/ less than 4 particles
    t.center_halos()
    t.get_covariance_matrices()
    t.get_eigenvectors()
    t.convert_bases()
    t.get_radii()    # center_halo(), get_covariance_matrices() and get_eigenvectors() functions must be called before
    t.get_half_mass_radii()
    return t        # t.halos is a list of halos contained inside the bgc2 file


def ascii_test(path='..\data\ellipsoid.dat'):
    """
    read data from an ascii file containing a single halo. By default, columns 1,2,3 (0-indexed) contain x,y,z coordinates,
    and first line of ascii file is skipped. See halos/helpers/helpers.py -> read_ascii_pos() for more documentation.
    """
    print "reading file: ", path
    s_time = time.clock()
    coords = helpers.read_ascii_pos(path)
    h = Halo('test', (0, 0, 0), coords)
    print 'execution time: ', time.clock()-s_time
    return h

# Generating a sample halo without ascii or bgc2 file
# x = [1,2,3,4,5]
# y = [2,3,5,2,1]
# z = [2,3,4,4,0]
# id = 'test_halo'
# center = (0, 0, 1)
# halo = helpers.create_halo(id, center, x, y, z)

# Perform calculations:
# halo.center_halo()
# halo.get_covariance_matrix()

# Or just pass it to do_all():
# halo = do_all(halo)

# Visualize halo:
# halo.visualize(ellipsoids=True)


def multicore_test(n):
    M = m.Multicore(n)
    M.get_data_from_files('../data/halos_0.1.bgc2')
    M.balance_load()
    #M.visualize()
    M.begin()
    r = M.get_results()
    print r
    print sum(r)

def axes_test():
    d1 = (1,3,4,5,3,2)
    d2 = (0.1, 0.5, 0.9, 0.6, 0.8, 0.4)
    bar_width = 0.35
    index = np.arange(len(d1))
    fig, ax1 = plt.subplots()
    ax1.bar(index, d1, bar_width, label='d1', color='b')
    ax1.set_ylabel('D1', color='b')
    for tl in ax1.get_yticklabels():
        tl.set_color('b')
    ax1.set_xticks(index+bar_width)
    ax1.set_xticklabels(('A', 'B', 'C', 'D', 'E', 'F'))
    ax1.set_title('Title')
    ax1.set_xlabel('X axis')

    ax2 = ax1.twinx()
    ax2.bar(index+bar_width, d2, bar_width, label='d2', color='r')
    ax2.set_ylabel('D2', color='r')
    for tl in ax2.get_yticklabels():
        tl.set_color('r')


    fig.show()

def worker(i, dataq, lock):
    interval = random.randint(1,5)
    lock.acquire()
    print 'Process ', os.getpid(), ' sleeping for ', interval, 's'
    lock.release()
    time.sleep(interval)
    dataq.put([i, os.getpid(), interval])


def master(n):
    L = mp.Lock()
    processes = []
    dataq = mp.Queue()
    for i in range(n):
        p = mp.Process(target=worker, args=(i, dataq, L))
        processes.append(p)
        p.start()
    print 'all processes started, getting data'
    for i in range(n):
        print dataq.get()
    for p in processes:
        p.join()
