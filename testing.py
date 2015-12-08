__author__ = 'Ibrahim'

from halos import *
import time


def main(path='..\data\halos_0.1.bgc2'):
    """
    reading data from a sample bgc2 file
    :param path: path to file
    :return: a HalfMassRadius instance containing list of halos
    """
    t = HalfMassRadius(path)
    t.read_data()
    t.filter(4)
    t.center_halos()
    t.get_covariance_matrices()
    t.get_eigenvectors()
    t.get_radii()
    t.get_half_mass_radii()
    return t


def ascii_test(path='..\data\ellipsoid.dat'):
    print "reading file: ", path
    s_time = time.clock()
    coords = read_ascii_pos(path)
    h = Halo('test', (0, 0, 0), coords)
    print 'execution time: ', time.clock()-s_time
    return h

def do_all(halo):
    s_time = time.clock()
    print 'beginning processing of halo:', halo.id
    halo.center_halo()
    halo.get_covariance_matrix()
    halo.get_eigenvectors()
    halo.get_radii()
    halo.get_half_mass_radius()
    print 'finished, execution time: ', time.clock()-s_time
    return halo
