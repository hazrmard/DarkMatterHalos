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
        for halo in self.halos:
            halo.cov = np.cov([halo.particles.x, halo.particles.y, halo.particles.z])
        print "covariance matrices obtained"

    def get_eigenvectors(self):
        for halo in self.halos:
            _, halo.eig = np.linalg.eig(halo.cov)


class Halo:
    coord_type = np.dtype([('x', np.float32), ('y', np.float32), ('z', np.float32)])

    def __init__(self, i, pos, particles):
        self.id = i
        self.pos = np.array(pos, dtype=Halo.coord_type).view(np.recarray)
        self.particles = particles
        self.cov = np.empty((3,3), dtype=np.float32)
        self.eig = np.empty((3,3), dtype=np.float32)
