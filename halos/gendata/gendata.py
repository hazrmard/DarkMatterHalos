__author__ = 'Ibrahim'

import numpy as np
from ..helpers import *


def _gen_angles(n=100):
    theta = np.random.uniform(0, np.pi, n)
    phi = np.random.uniform(0, 2 * np.pi, n)
    return theta, phi


def _redestribute(a, pdf=lambda x: x):
    '''Redestributes range of points according to the inverse of the  probability distribution function/
    :a array of points (like random uniform samples)
    :pdf function that takes 'a' and returns another list. e.g. def (x): return [i**2 for  i in x]
    '''
    min_a = min(a)
    max_a = max(a)
    range_a = max_a - min_a
    b = pdf(a)
    min_b = min(b)
    max_b = max(b)
    range_b = max_b - min_b
    b_centered = b - min_b
    b_normed = b_centered / range_b
    a_new = b_normed * range_a
    a_new = min_a + a_new
    return a_new


def shell(n=100, r=(10, 10, 10)):
    theta, phi = _gen_angles(n)
    x = r[0] * np.sin(theta) * np.cos(phi)
    y = r[1] * np.sin(theta) * np.sin(phi)
    z = r[2] * np.cos(theta)
    return x, y, z


def two_shells(n=(100, 100), r=(10, 10, 10), scale=0.5):
    x, y, z = shell(n[0], (r[0]*scale, r[1]*scale, r[2]*scale))
    x2, y2, z2 = shell(n[1], r)
    x = np.append(x, x2)
    y = np.append(y, y2)
    z = np.append(z, z2)
    return x, y, z

    
def random_rect(n=100, dims=(10,10,10), center=(0,0,0)):
    x = np.random.uniform(-dims[0], dims[0], (1,n))[0] - center[0]
    y = np.random.uniform(-dims[1], dims[1], (1,n))[0] - center[1]
    z = np.random.uniform(-dims[2], dims[2], (1,n))[0] - center[2]
    id = 'random_rect'
    return helpers.create_halo(id, center, x, y, z)

def random_halo(n=100, dims=(10,10,10), center=(0,0,0), pdf=lambda x: [i**2 for i in x]):
    '''Generate a halo of particles following some particle distribution.
    :n number of particles
    :dims tuple of radial size in each direction
    :center location of center of halo
    :pdf probability density function for particle distribution. Should accept and return a list of numbers.
    '''
    v, u = _gen_angles(n=n)
    rx = _redestribute(np.random.uniform(0, dims[0], n), pdf= pdf)
    ry = _redestribute(np.random.uniform(0, dims[1], n), pdf= pdf)
    rz = _redestribute(np.random.uniform(0, dims[2], n), pdf= pdf)
    x = rx * np.cos(u) * np.sin(v) + center[0]
    y = ry * np.sin(u) * np.sin(v) + center[1]
    z = rz * np.cos(v) + center[2]
    id = 'random_halo'
    return helpers.create_halo(id, center, x, y, z)
    