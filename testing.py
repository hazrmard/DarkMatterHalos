__author__ = 'Ibrahim'

# This file contains examples of functions that can be used to process halo data.

from halos import *
import time


def bgc2_test(path='..\data\halos_0.1.bgc2'):
    """
    reading data from a sample bgc2 file containing multiple halos
    :param path: path to file
    :return: a HalfMassRadius instance containing list of halos (<return_variable>.halos)
    """
    t = HalfMassRadius(path)
    t.read_data()
    t.filter(4)         # filter out halos w/ less than 4 particles
    t.center_halos()
    t.get_covariance_matrices()
    t.get_eigenvectors()
    t.get_radii()    # center_halo(), get_covariance_matrices() and get_eigenvectors() functions must be called before
    t.get_half_mass_radii()
    return t


def ascii_test(path='..\data\ellipsoid.dat'):
    """
    read data from an ascii file containing a single halo. By default, columns 1,2,3 (0-indexed) contain x,y,z coordinates,
    and first line of ascii file is skipped. See halos/helpers/helpers.py -> read_ascii_pos() for more documentation.
    """
    print "reading file: ", path
    s_time = time.clock()
    coords = read_ascii_pos(path)
    h = Halo('test', (0, 0, 0), coords)
    print 'execution time: ', time.clock()-s_time
    return h

def do_all(halo):
    """
    take a halo instance and take it through all functions necessary for calculating half mass radius.
    Returns halo instance with finished calculations.
    """
    s_time = time.clock()
    print 'beginning processing of halo:', halo.id
    halo.center_halo()
    halo.get_covariance_matrix()
    halo.get_eigenvectors()
    halo.get_radii()
    halo.get_half_mass_radius()
    print 'finished, execution time: ', time.clock()-s_time
    return halo

# Generating a sample halo without ascii or bgc2 file
# x = [1,2,3,4,5]
# y = [2,3,5,2,1]
# z = [2,3,4,4,0]
# id = 'test_halo'
# center = (0, 0, 1)
# halo = create_halo(id, center, x, y, z)

# Perform calculations:
# halo.center_halo()
# halo.get_covariance_matrix()

# Or just pass it to do_all():
# halo = do_all(halo)

# Visualize halo:
# halo.visualize(ellipsoids=True)
