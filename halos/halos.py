#!/usr/bin/env python

import sys
import numpy as np
import warnings
import matplotlib.pyplot as plt
import bgc2
import config
from mpl_toolkits.mplot3d import Axes3D


#Settings
#plt.hold(True)
warnings.simplefilter('error', RuntimeWarning)      # raise exception on RuntimeWarning

COORDS = config.COORDS

class HalfMassRadius:
    def __init__(self, file):
        self.file = file
        self.header = []
        self.halos = []     # list of halos
        self.h = []
        self.warnings = {}

    def read_data(self):
        """
        bcg2.read_bcg2_numpy returns 3 numpy Record Arrays for header, halos, and particles.
        Header and halos are single dimensional Record Arrays containing data and halo information.
        Particle is a Record Array of Record Arrays for each halo id and for each particle axis.
        Schema for arrays can be found in halos/helpers/bcg2.py.
        """
        self.header, halos, particles = bgc2.read_bgc2_numpy(self.file)
        self.h = halos
        for i in xrange(len(halos)):
            self.halos.append(Halo(halos[i].id, (halos[i].x, halos[i].y, halos[i].z), particles[i]))
        print "data file read"

    def filter(self, minimum=None, maximum=None):
        if minimum is not None:
            self.halos = [h for h in self.halos if len(h.particles) >= minimum]
            print "halos with less than " + str(minimum) + " particles filtered out"
        if maximum is not None:
            self.halos = [h for h in self.halos if len(h.particles) <= maximum]
            print "halos with more than " + str(maximum) + " particles filtered out"

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
    
    def convert_bases(self):
        for halo in self.halos:
            halo.convert_basis()
        print "bases converted"

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
    coord_type = COORDS 

    def __init__(self, i, pos, particles):
        self.id = i
        self.pos = np.array(pos, dtype=Halo.coord_type).view(np.recarray)
        self.particles = particles  # record array of particles (see halos/helpers/helpers.py -> create_halo())
        self.particlesn = particles     # record array of particles transformed to new basis
        self.cov = np.empty((3, 3), dtype=np.float32)   # covariance matrix
        self.eig = np.empty((3, 3), dtype=np.float32)    # eigenvectors
        self.evals = np.empty(3, dtype=np.float32)       # eigenvalues
        self.radii = np.empty(len(self.particles), dtype=np.float32)    # ellipsoidal radii of particles
        self.half_mass_radius = np.float32(0)                                  # ellipsoidal radius of hald mass sub-halo
        self.fig = None

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
        self.eig = np.vstack((eigvecs[:, order[0]], eigvecs[:, order[1]], eigvecs[:, order[2]])).T
    
    def convert_basis(self, basis=None):
        """
        Transform particles to new basis defined by eigenvectors.
        :basis 3 x 3 np array of new basis column vectors
        """
        basis = self.eig if basis is None else basis
        coords = np.array(zip(self.particles.x, self.particles.y, self.particles.z))    # n x 3 matrix
        e_1 = np.linalg.inv(basis)                     # inverse eigenvector matrix
        e_1c = np.dot(e_1, coords.T)               # points transformed to the new basis.
        self.particlesn = np.array(zip(*e_1c), dtype=COORDS).view(np.recarray)

    def get_radii_2(self):          # depracated approach to finding ellipsoid fit. Use get_radii()
        """
        Using the general ellipsoid equation x(T).A.x = er^2 (ellipsoidal radius) to find 'er' for all particles.
        From Wikipedia: the eigenvectors of 'A' define the ellipsoid's principal components, therefore A is either
        the covariance matrix or its inverse. My calculations on Mathematica indicate the latter.
        """
        coords = np.array(zip(self.particles.x, self.particles.y, self.particles.z))    # n x 3 matrix
        invcov = np.linalg.inv(self.cov)
        invdet = np.linalg.det(self.cov)
        self.radii = np.sqrt(np.einsum('ij,ji->i', coords, np.dot(invcov / invdet, coords.T)))

    def get_radii(self, evals=None):
        """
        The points are first transformed into a basis defined by the eigenvectors of the covariance matrix.
        Using equation ax^2 + by^2 + cz^2 = er where a, b, and c are normalized eigenvalues of the covariance matrix.
        """
        evals = self.evals if evals is None else evals
        coords = np.array(zip(self.particlesn.x, self.particlesn.y, self.particlesn.z)).T    # 3 x n matrix
        c2 = coords * coords                        # squared coordinates. In column vectors ([x^2],[y^2],[z^2])
        w = np.diag(evals) / max(evals)        # normalized diagonal matrix containing eigenvalues of covariance matrix
        wc2 = np.dot(w, c2)                          # weighted squared coords. In column vectors ([ax^2],[by^2],[cz^2])
        self.radii = np.sqrt(np.einsum('ij->j', wc2))    # sqrt of sum of weighted squared coordinates

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
    
    def higher_order_fit(self, order=2):
        """
        Create sub-halo from particles inside half-mass radius, recompute covariance matrix, recalculate particles
        inside half-mass radius, repeat.
        :order number of times -1 to repeat
        """
        for i in range(1, order):
            indices, _ = self.cut()
            coords =  zip(self.particlesn.x[indices], self.particlesn.y[indices], self.particlesn.z[indices])
            h = Halo(i='subhalo', pos=(0,0,0), particles=np.array(coords, dtype=self.coord_type).view(np.recarray))
            h.get_covariance_matrix()
            h.get_eigenvectors()
            self.get_radii(evals=h.evals)
    
    def encapsulation(self):
        """
        Get a percentage of points that are not encapsulated by the ellipsoid, but should be.
        """
        indices, _ = self.cut()
        r = self.cleave(indices)
        norm_evals = self.evals / np.amax(self.evals)
        R = np.sqrt(np.sum(np.square(norm_evals * np.array([r.x, r.y, r.z]))))
        total = np.sum(np.where(self.radii[indices]>R, 1, 0))
        return total / float(len(indices))
    
    def cut(self, fraction=1.0):
        """
        Returns indices of the first and second halves of particle/radius arrays corresponding to particles
        in and outside the half mass radius.
        :fraction fraction of particles to select from each cut
        """
        sortedindices = np.array(np.argsort(self.radii))
        l = int(len(sortedindices)/2)
        firsthalf = np.random.choice(sortedindices[:l], size=int(l*fraction))
        secondhalf = np.random.choice(sortedindices[l:], size=len(sortedindices)-int(l*fraction))
        return firsthalf, secondhalf
    
    def cleave(self, indices=None):
        """
        Returns radii along each principal component corresponding the largest projection on that component
        by particles in the list of indices. The the fitted ellipsoids will cleave to the outermost particle.
        :indices List of indices for which to get the radii.
        returns a numpy record array of radii e.g. return.x, return.y, return.z give each radius
        """
        indices = list(range(len(self.radii))) if indices is None else indices
        Rx = np.amax(np.absolute(self.particlesn.x[indices]))
        Ry = np.amax(np.absolute(self.particlesn.y[indices]))
        Rz = np.amax(np.absolute(self.particlesn.z[indices]))
        return np.array((Rx, Ry, Rz), dtype=Halo.coord_type).view(np.recarray)

    def visualize(self, ellipsoids=False, fraction=1.0):
        """
        3D plot of particles. Particles within half mass radius are in red. Others are in blue.
        :fraction fraction of particles to show (1.0 means 100%)
        :ellipsoids whether to draw fitted surfaces
        """
        self.fig = plt.figure(self.id).add_subplot(111, projection='3d')
        firsthalf, secondhalf = self.cut(fraction)
        self.fig.scatter(self.particlesn.x[firsthalf], self.particlesn.y[firsthalf], self.particlesn.z[firsthalf], c='r')
        self.fig.scatter(self.particlesn.x[secondhalf], self.particlesn.y[secondhalf], self.particlesn.z[secondhalf], c='b')
        if ellipsoids:
            self._draw_ellipsoids((firsthalf, secondhalf))
        plt.show()
        plt.close()

    def _draw_ellipsoids(self, indices):
        ax = self.fig
        radii = [self.cleave(indices[0]), self.cleave(np.append(indices[1], indices[0]))]
        colors = ('r', 'c')
        alpha = 0.3
        
        u = np.linspace(0, 2 * np.pi, 100)              # generate ellipsoid surface data parametrically
        v = np.linspace(0, np.pi, 100)
        
        for r, c in zip(radii, colors):            
            x = r.x * np.outer(np.cos(u), np.sin(v))
            y = r.y * np.outer(np.sin(u), np.sin(v))
            z = r.z * np.outer(np.ones_like(u), np.cos(v))
            tcoords = np.dot(np.identity(3), np.array(zip(x, y, z)).T)        # ellipsoid coordinates
            ax.plot_surface(tcoords[0, :], tcoords[1, :], tcoords[2, :], rstride=4, cstride=4, color=c, alpha=alpha)
            alpha /= 2.0
        max_r = np.amax([radii[1].x, radii[1].y, radii[1].z])
        ax.set_xlim3d(-max_r, max_r)   # set axes limits
        ax.set_ylim3d(-max_r, max_r)
        ax.set_zlim3d(-max_r, max_r)

