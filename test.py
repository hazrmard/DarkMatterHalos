__author__ = 'Ibrahim'

from halos import *

filepath = ""
test = HalfMassRadius(filepath)
test.read_data()
test.center_halos()
test.get_covariance_matrices()
test.get_radii()
test.get_half_mass_radii()