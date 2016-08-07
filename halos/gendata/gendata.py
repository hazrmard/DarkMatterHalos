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


def shell(n=100, dims=(10, 10, 10), center=(0,0,0)):
    theta, phi = _gen_angles(n)
    x = dims[0] * np.sin(theta) * np.cos(phi)
    y = dims[1] * np.sin(theta) * np.sin(phi)
    z = dims[2] * np.cos(theta)
    return helpers.create_halo('shell', (0,0,0), x, y, z)


def two_shells(n=(100, 100), dims=(10, 10, 10), center=(0,0,0), scale=0.5):
    h1 = shell(n=n[0], center=(dims[0]*scale, dims[1]*scale, dims[2]*scale))
    h2 = shell(n[1], (dims[0], dims[1], dims[2]))
    x = np.append(h1.particles.x, h2.particles.x) + center[0]
    y = np.append(h1.particles.y, h2.particles.y) + center[1]
    z = np.append(h1.particles.z, h2.particles.z) + center[2]
    return helpers.create_halo('two-shells', (0,0,0), x, y, z)


def rect(n=100, dims=(10,10,10), center=(0,0,0)):
    x = np.random.uniform(-dims[0], dims[0], (1,n))[0] - center[0]
    y = np.random.uniform(-dims[1], dims[1], (1,n))[0] - center[1]
    z = np.random.uniform(-dims[2], dims[2], (1,n))[0] - center[2]
    id = 'random_rect'
    return helpers.create_halo(id, center, x, y, z)


def halo(n=100, dims=(10,10,10), center=(0,0,0), pdf=lambda x: [np.sqrt(i) for i in x]):
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
