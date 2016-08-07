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
    :return: a Halos instance containing list of halos (<return_variable>.halos)
    """
    t = Halos(path)
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
    do_all(h)
    return h

# Generating a sample halo without ascii or bgc2 file
# ids = [0,1,2,3,4]
# x = [1,2,3,4,5]
# y = [2,3,5,2,1]
# z = [2,3,4,4,0]
# id = 'test_halo'
# center = (0, 0, 1)
# halo = helpers.create_halo(id, center, x, y, z, ids)

# Perform calculations:
# halo.center_halo()
# halo.get_covariance_matrix()

# Or just pass it to do_all():
# halo = do_all(halo)

# Visualize halo:
# halo.visualize(ellipsoids=True)

if __name__=='__main__':
    print('Running calculations on sample BGC2 file:')
    bgc2_test()
    print('\nRunning calculations on sample ascii file:')
    ascii_test()
