__author__ = 'Ibrahim'

import numpy as np


def _gen_angles(n=100):
    theta = np.arccos(np.random.uniform(-1, 1, n))
    phi = np.random.uniform(0, 2 * np.pi, n)
    return theta, phi


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
