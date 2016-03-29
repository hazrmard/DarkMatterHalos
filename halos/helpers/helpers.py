import numpy as np
from .. import config
from .. import halos


def read_ascii_pos(filepath='', settings=config._default_ascii_settings):
    """
    this function reads columns from an ascii file into a numpy array
    :param filepath: path of file
    :param settings: a dictionary of arguments to be applied to the file. Argument documentation can be seen @
    http://docs.scipy.org/doc/numpy/reference/generated/numpy.genfromtxt.html
    :return: a numpy record array
    """
    data = np.genfromtxt(filepath, **settings)
    data = data.view(np.recarray)
    return data


def create_halo(haloid, halocenter, halox, haloy, haloz):
    """
    given array of x,y,z coordinates generates a Halo instance
    :param haloid: id of halo
    :param halocenter: centre of halo as a touple e.g (0,0,0)
    :param halox: list/array of x coordinates
    :param haloy: list/array of y coordinates
    :param haloz: list/array of z coordinates
    :return: a Halo instance
    """
    coords = zip(halox, haloy, haloz)
    h = halos.Halo(haloid, halocenter, np.array(coords, dtype=config.COORDS).view(np.recarray))
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
    halo.convert_basis()
    halo.get_radii()
    halo.get_half_mass_radius()
    print 'finished, execution time: ', time.clock()-s_time
    return halo

def std_dev(arr):
    n = len(arr)
    ssum = sum([x**2 for x in arr])
    return ssum/n - (sum(arr)/n)**2
