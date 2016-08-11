__author__ = 'Ibrahim'

# This file contains examples of functions that can be used to process halo data.

from halos import Halos, gendata, helpers


def bgc2_test(path='data\\halos_0.1.bgc2'):
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

def bgc2_merger_test(f1=r'data\*0000.bgc2', f2=r'data\*0001.bgc2'):
    try:
        h = Halos(f1, verbose=False)
        g = Halos(f2, verbose=False)
        h.read_data(level=0)
        g.read_data(level=0)
    except Exception as e:
        print('Error reading test data files: data\\*0000.bgc2 and data\\*0001.bgc2')
        print(str(e))
        return
    try:
        h+g
        print('SUCCESS: Compatible files merged.')
    except:
        print('FAILURE: BGC2 merger failed.')
    try:
        h+h
    except ValueError as e:
        print('SUCCESS: Invalid merger successfully prevented.')

def ascii_test(path='data\\ellipsoid.dat'):
    """
    read data from an ascii file containing a single halo. By default, columns 1,2,3 (0-indexed) contain x,y,z coordinates,
    and first line of ascii file is skipped. See halos/helpers/helpers.py -> read_ascii_pos() for more documentation.
    """
    print "reading file: ", path
    coords = helpers.read_ascii_pos(path)
    h = Halo('test', (0, 0, 0), coords)
    helpers.do_all(h)
    return h

def random_halo_test():
    h = gendata.halo(n=1000)        # generate random halo
    helpers.do_all(h)               # all the core functions to find half-mass r
    h.higher_order_fit(order=5)
    h.report()
    h.visualize(ellipsoids=True)

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
# halo = helpers.do_all(halo)

# Visualize halo:
# halo.visualize(ellipsoids=True)

if __name__=='__main__':
    print('Running calculations on random halos:')
    random_halo_test()

    try:
        print('\nRunning calculations on sample BGC2 file:')
        bgc2_test()
    except IOError as e:
            print('Test bgc2 file \'data\\halos_0.1.bgc2\' not found.')

    print('Testing BGC2 file merger with Halos class:')
    bgc2_merger_test()

    try:
        print('\nRunning calculations on sample ascii file:')
        ascii_test()
    except IOError as e:
            print('Test ascii file \'data\\ellipsoids.dat\' not found.')
