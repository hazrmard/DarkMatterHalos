#!/usr/bin/env python

import sys
import bgc2
import numpy as np

class HalfMassRadius:
    def __init__(self, file):
        self.file = file
        self.header = None
        self.halos = None
        self.particles = None

    def read_data(self):
        """
            bcg2.read_bcg2_numpy returns 3 numpy Record Arrays for header, halos, and particles.
            Header and halos are single dimensional Record Arrays containing data and halo information.
            Particle is a Record Array of Record Arrays for each halo id and for each particle.
            Schema for arrays can be found in bcg2.py.
        """
        self.header, self.halos, self.particles = bgc2.read_bgc2_numpy(self.file)
        print "data file read"

    def center_halos(self):
        """
            Iterating over each halo to obtain its centre. Then iterating over each particle in the halo to center it.
        """
        for i in xrange(len(self.halos)):       # iterating over each halo
            hx = self.halos[i].x
            hy = self.halos[i].y
            hz = self.halos[i].z
            for j in xrange(len(self.particles[i])):    # iterating over all particles/halo and centering them
                self.particles[i][j].x -= hx
                self.particles[i][j].y -= hy
                self.particles[i][j].z -= hz

