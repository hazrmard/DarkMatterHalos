#!/usr/bin/env python

import sys
import bgc2
import numpy as np


class HalfMassRadius:
    def __init__(self, file):
        self.file = file
        self.header = []
        self.halos = []
        self.h = []

    def read_data(self):
        """
        bcg2.read_bcg2_numpy returns 3 numpy Record Arrays for header, halos, and particles.
        Header and halos are single dimensional Record Arrays containing data and halo information.
        Particle is a Record Array of Record Arrays for each halo id and for each particle.
        Schema for arrays can be found in bcg2.py.
        """
        self.header, halos, particles = bgc2.read_bgc2_numpy(self.file)
        self.h = halos
        for i in xrange(len(halos)):
            self.halos.append(Halo(halos[i].id, (halos[i].x, halos[i].y, halos[i].z), particles[i]))
        print "data file read"

    def center_halos(self):
        """
        Iterating over each halo and subtracting its centre position from particle coordinates.
        """
        for halo in self.halos:
            halo.particles.x -= halo.pos.x
            halo.particles.y -= halo.pos.y
            halo.particles.z -= halo.pos.z
        print "halos centered"

    def get_covariance_matrices(self):
        """
        Calculate the covariance matrix for each halo particle coordinate set.
        """
        for halo in self.halos:
            halo.cov = np.cov([halo.particles.x, halo.particles.y, halo.particles.z])
        print "covariance matrices obtained"

    def get_eigenvectors(self):
        """
        Calculate eigenvectors for each halo and store them in order of decreasing eigenvalue. The first Principal
        Component (major axis) is in the first column.
        """
        for halo in self.halos:
            eigenvals, eigvecs = np.linalg.eig(halo.cov)
            order = np.argsort(eigenvals)[::-1]
            halo.eig = np.vstack((eigvecs[:order[0]], eigvecs[:order[1]], eigvecs[:order[2]])).T

    def get_radii(self):
        """
        Using the general ellipsoid equation x(T).A.x = ellipsoidal radius (er) to find 'er' for all particles.
        From Wikipedia: the eigenvectors of 'A' define the ellipsoid's principal components, therefore A is the
        covariance matrix.
        """
        for halo in self.halos:
            for i in xrange(len(halo.particles)):
                coords = np.vstack((halo.particles[i].x, halo.particles[i].y, halo.particles[i].z))
                er = np.dot(coords.T, np.dot(halo.cov, coords))
                np.append(halo.radii, er)
        print "radii computed"

    def get_half_mass_radii(self):
        """
        Sorting radii by size and taking the median of the array as half_mass_radius.
        """
        for halo in self.halos:
            radii = np.sort(halo.radii)
            l = len(radii)
            if len(halo.radii) % 2 == 0:
                halo.half_mass_radius = np.float32((radii[l/2] + radii[(l-2)/2])/2.0)
            else:
                halo.half_mass_radius = np.float32(radii[(l-1)/2.0])
        print "half mass radii computed"


class Halo:
    coord_type = np.dtype([('x', np.float32), ('y', np.float32), ('z', np.float32)])

    def __init__(self, i, pos, particles):
        self.id = i
        self.pos = np.array(pos, dtype=Halo.coord_type).view(np.recarray)
        self.particles = particles
        self.cov = np.empty((3, 3), dtype=np.float32)
        self.eig = np.empty((3, 3), dtype=np.float32)
        self.radii = np.empty(len(self.particles), dtype=np.float32)
        self.half_mass_radius = np.float32(0)
