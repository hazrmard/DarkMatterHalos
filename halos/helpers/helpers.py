import numpy as np
from .. import config
from .. import halos
import glob
import os
import re
from operator import itemgetter


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


def create_halo(haloid, halocenter, halox, haloy, haloz, partids=None):
    """
    given array of x,y,z coordinates generates a Halo instance
    :param haloid: id of halo
    :param halocenter: centre of halo as a touple e.g (0,0,0)
    :param partids: list/array of integer id values for each particle
    :param halox: list/array of x coordinates
    :param haloy: list/array of y coordinates
    :param haloz: list/array of z coordinates
    :return: a Halo instance
    """
    if partids is None:
        partids = range(len(halox))
    coords = zip(partids, halox, haloy, haloz)
    h = halos.Halo(haloid, halocenter, np.array(coords, dtype=config.COORDSwID).view(np.recarray))
    return h

def do_all(halo):
    """
    take a halo instance and take it through all functions necessary for calculating half mass radius.
    Returns halo instance with finished calculations.
    """
    #s_time = time.clock()
    #print 'beginning processing of halo:', halo.id
    halo.center_halo()
    halo.get_covariance_matrix()
    halo.get_eigenvectors()
    halo.convert_basis()
    halo.get_radii()
    halo.get_half_mass_radius()
    #print 'finished, execution time: ', time.clock()-s_time
    return halo

def std_dev(arr):
    n = len(arr)
    ssum = sum([x**2 for x in arr])
    return ssum/n - (sum(arr)/n)**2

def generate_file_groups(path, version_levels=1, ignore_ext=True):
    """create a list of lists, where each inner list is a group of files with
    the same version level in the file name. Assumes file naming conventions are
    consistent. Intended for BGC2 snapshots which are stored over multiple files.
    :param path: UNIX style path i.e. *,?,[] wildcards accepted
    :param version_levels: which subversion level to group on. For e.g filenames:
    F1.0.0, F1.0.1, F1.1,0, F1.1.1, with version_levels=0 each file will have its
    own group. With version_levels=1, [F1.0.0, F1.0.1], [F1.1,0, F1.1.1] => 2
    groups w/ vers# 1.0.* and 1.1.*. And with version_levels=2, 1 group of vers# 1.*.*
    :param ignore_ext: Whether or not to include any numbers in extension i.e
    BGC'2' which is not a version indicator. True removes file extension.
    """
    flist = glob.glob(path)
    finfo = []
    if len(flist)==0:
        raise Exception("No files found.")

    if ignore_ext:      # remove extension (chars inc. and after last .)
        finfo = [os.path.splitext(f)[0] for f in flist] # filepath w/0 extension
    else:
        finfo = [x for x in flist]

    pattern = re.compile(r'\d+')
    finfo = [[int(num) for num in pattern.findall(f)] for f in finfo] # nums in filepath representing ver#
    possible_subversions = len(finfo[0])
    if (possible_subversions-1<version_levels):
        raise ValueError('Not enough subversions in filename for level ' + \
                            str(version_levels) + ' grouping.')

    names_and_info = zip(flist, finfo)  # file paths combined with extracted ver #s
    sort_indices = range(0, possible_subversions-1)    # ver # indices to sort by
    key = itemgetter(*sort_indices)    # to get subversions from ver #
    names_and_info.sort(key=lambda x: key(x[1]))

    groups=[]
    ngroups=0
    i=0
    while i< len(names_and_info):
        groups.append([])
        curr_ver = names_and_info[i][1][0:possible_subversions-version_levels]
        while (i<len(names_and_info) and \
            curr_ver==names_and_info[i][1][0:possible_subversions-version_levels]):
            groups[ngroups].append(names_and_info[i][0])
            i+=1
        ngroups+=1

    return groups
