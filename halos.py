#!/usr/bin/env python

import sys
import bgc2
import numpy as np
import warnings
warnings.simplefilter('error', RuntimeWarning)


class HalfMassRadius:
    def __init__(self, file):
        self.file = file
        self.header = []
        self.halos = []
        self.h = []
        self.warnings = {}

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
        for halo in self.halos:
            halo.center_halo()
        print "halos centered"

    def get_covariance_matrices(self):
        for halo in self.halos:
            halo.get_covariance_matrix()
        print "covariance matrices obtained"

    def get_eigenvectors(self):
        for halo in self.halos:
            halo.get_eigenvectors()
        print "eigenvectors computed"

    def get_radii(self):
        for halo in self.halos:
            try:
                halo.get_radii()
            except np.linalg.linalg.LinAlgError:        # singular covariance matrix
                halo.radii = -1                         # handle singular covariance matrices
            except RuntimeWarning as a:
                self.warnings[halo.id] = a.message
        print "total warnings: " + str(len(self.warnings))
        print "radii computed"

    def get_half_mass_radii(self):
        for halo in self.halos:
            halo.get_half_mass_radius()
        print "half mass radii computed"


class Halo:
    coord_type = np.dtype([('x', np.float32), ('y', np.float32), ('z', np.float32)])

    def __init__(self, i, pos, particles):
        self.id = i
        self.pos = np.array(pos, dtype=Halo.coord_type).view(np.recarray)
        self.particles = particles
        self.cov = np.empty((3, 3), dtype=np.float32)
        self.eig = np.empty((3, 3), dtype=np.float32)
        self.evals = np.empty(3, dtype=np.float32)
        self.radii = np.empty(len(self.particles), dtype=np.float32)
        self.half_mass_radius = np.float32(0)

    def center_halo(self):
        """
        Iterating over each halo and subtracting its centre position from particle coordinates.
        """
        self.particles.x -= self.pos.x
        self.particles.y -= self.pos.y
        self.particles.z -= self.pos.z

    def get_covariance_matrix(self):
        """
        Calculate the covariance matrix for each halo particle coordinate set.
        """
        self.cov = np.cov([self.particles.x, self.particles.y, self.particles.z])

    def get_eigenvectors(self):
        """
        Calculate eigenvectors for each halo and store them in order of decreasing eigenvalue. The first Principal
        Component (major axis) is in the first column.
        """
        eigenvals, eigvecs = np.linalg.eig(self.cov)
        order = np.argsort(eigenvals)[::-1]
        self.evals = np.sort(eigenvals)[::-1]
        self.eig = np.vstack((eigvecs[:order[0]], eigvecs[:order[1]], eigvecs[:order[2]])).T

    def get_radii(self):
        """
        Using the general ellipsoid equation x(T).A.x = er^2 (ellipsoidal radius) to find 'er' for all particles.
        From Wikipedia: the eigenvectors of 'A' define the ellipsoid's principal components, therefore A is either
        the covariance matrix or its inverse. My calculations on Mathematica indicate the latter.
        """
        coords = np.array(zip(self.particles.x, self.particles.y, self.particles.z))    # n x 3 matrix
        self.radii = np.sqrt(np.einsum('ij,ji->i', coords, np.dot(np.linalg.inv(self.cov), coords.T)))

    def get_half_mass_radius(self):
        """
        Sorting radii by size and taking the median of the array as half_mass_radius.
        """
        radii = np.sort(self.radii)
        l = len(radii)
        if len(self.radii) % 2 == 0:
            self.half_mass_radius = np.float32((radii[int(l/2)] + radii[(l-2)/2])/2.0)
        else:
            self.half_mass_radius = np.float32(radii[int((l-1)/2.0)])
