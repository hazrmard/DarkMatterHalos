import numpy as np
from ..halos import Halo

COORDS = np.dtype([('x', np.float32), ('y', np.float32), ('z', np.float32)])    # coordinates data type

_default_ascii_settings = {'dtype':COORDS,
                           'comments':'#',
                           'delimiter':None,
                           'skip_header': 1,
                           'skip_footer':0,   #   cannot be used together with 'max_rows' argument
                           'converters':None,
                           'missing_values':None,
                           'filling_values':None,
                           'usecols': (1, 2, 3),
                           'names':None,
                           'excludelist':None,
                           'deletechars':None,
                           'replace_space':'_',
                           'autostrip':False,
                           'case_sensitive':True,
                           'defaultfmt':'f%i',
                           'unpack':None,
                           'usemask':False,
                           'loose':True,
                           'invalid_raise':True,
                           # 'max_rows':200
                           }


def read_ascii_pos(filepath='', settings=_default_ascii_settings):
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
    h = Halo(haloid, halocenter, np.array(coords, dtype=COORDS).view(np.recarray))
    return h
