#!/usr/bin/env python

import sys
import numpy as np
import config
import warnings
import matplotlib.pyplot as plt
import bgc2
from mpl_toolkits.mplot3d import Axes3D
import glob
plt.ioff()

#Settings
#plt.hold(True)
warnings.simplefilter('error', RuntimeWarning)      # raise exception on RuntimeWarning

COORDS = config.COORDS

class Halos:
    def __init__(self, files, verbose=True):
        """
        :param files: single or list of paths. Paths can contain wildcards. Accepts
                      UNIX style path i.e. *,?,[] wildcards accepted
        """
        if type(files) is list:
            self.files = files
        else:
            self.files = [files]
        temp_files = []
        for file in self.files:
            temp_files.extend(glob.glob(file))
        self.files = temp_files
        if len(self.files)==0:
            raise Exception('No files found.')

        self.header = []
        self.halos = []     # list of halo objects w/ particle data
        self.h = []         # halo metadata (position, id, size etc.)
        self.warnings = {}
        self.verbose = verbose

    def read_data(self, level=2, sieve=None, strict=True, onlyid=False):
        """
        bcg2.read_bcg2_numpy returns 3 numpy Record Arrays for header, halos, and particles.
        Header and halos are single dimensional Record Arrays containing data and halo information.
        Particle is a Record Array of Record Arrays for each halo id and for each particle axis.
        Schema for arrays can be found in halos/helpers/bcg2.py.
        :level reads either header(0), halo metadata(1) or particle data (2). Should be 2.
        :sieve a set of ids to keep. All others discarded. If None, all data are kept.
        :strict[DEPRACATED]=False creates empty halo instances just with metadata, True does not create instances
        :onlyid=True reads only halo and particle ids, default is False for full data
        """
        for file in self.files:
            header, h, particles = bgc2.read_bgc2_numpy(file, level=level, sieve=sieve, onlyid=onlyid)

            if level>=0:                # appending header data, level=0
                self._add_header(header)

            if level>=1:                # appending halo group meta data, level=1
                self.h.extend(h)

            if level>=2:                # appending particle data, level=2
                if not onlyid:
                    for i in xrange(len(h)):
                        self.halos.append(Halo(h[i].id, (h[i].x, h[i].y, h[i].z), particles[i]))
                else:
                    for i in xrange(len(h)):
                        self.halos.append(Halo(h[i], (None, None, None), particles[i]))
            # if level==1 and strict==False:
            #     for i in xrange(len(h)):
            #         self.halos.append(Halo(h[i].id, (h[i].x, h[i].y, h[i].z), ()))
        if self.verbose:
            print "data file(s) read"

    def filter(self, minimum=None, maximum=None):
        """filter out halos with less/greater than # of particles
        :param minimum: halos with <= this # of particles will be removed
        :param maximim: halos with >= this # of particles will be removed
        """
        if minimum is not None:
            self.halos = [self.halos[i] for i in range(len(self.h)) if self.h[i].npart >= minimum]
            self.h = [h for h in self.h if h.npart >= minimum]
            if self.verbose:
                print "halos with less than " + str(minimum) + " particles filtered out"
        if maximum is not None:
            self.halos = [self.halos[i] for i in range(len(self.h)) if self.h[i].npart <= maximum]
            self.h = [h for h in self.h if h.npart <= maximum]
            if self.verbose:
                print "halos with more than " + str(maximum) + " particles filtered out"

    def center_halos(self):
        for halo in self.halos:
            halo.center_halo()
        if self.verbose:
            print "halos centered"

    def get_covariance_matrices(self):
        for halo in self.halos:
            halo.get_covariance_matrix()
        if self.verbose:
            print "covariance matrices obtained"

    def get_eigenvectors(self):
        for halo in self.halos:
            halo.get_eigenvectors()
        if self.verbose:
            print "eigenvectors computed"

    def convert_bases(self):
        for halo in self.halos:
            halo.convert_basis()
        if self.verbose:
            print "bases converted"

    def get_radii(self):
        for halo in self.halos:
            try:
                halo.get_radii()
            except np.linalg.linalg.LinAlgError:        # singular covariance matrix
                halo.radii = -1                         # handle singular covariance matrices
            except RuntimeWarning as a:
                self.warnings[halo.id] = a.message
        if self.verbose:
            print "total warnings: " + str(len(self.warnings))
            print "radii computed"

    def get_half_mass_radii(self):
        for halo in self.halos:
            _ = halo.get_half_mass_radius()
        if self.verbose:
            print "half mass radii computed"

    def higher_order_fit(self, order=2):
        i=0
        for halo in self.halos:
            try:
                halo.higher_order_fit(order=order)
            except:
                if self.verbose:
                    print "Complex warning with halo id:" + str(halo.id)
                i+=1
        if self.verbose:
            print "higher order fit completed, total warnings: " + str(i)

    def _add_header(self, header):
        """append other BGC2 file header and check for incompatibilities"""
        if len(self.header)==0:
            self.header = [header]
        elif self.header[0]==header:
            raise ValueError('Cannot merge. Multiple headers are identical.')
        else:
            # check for fields that must be equal
            for field in config.MANDATORY_HEADER_FIELDS:    # check compartibility
                if self.header[0][field]!=header[field]:
                    raise ValueError('Cannot merge. Field: ' + field + 'is not equal for both files.')
            # check for fields that must not be equal
            for i in range(len(self.header)):
                for field in config.EXCLUSIVE_HEADER_FIELDS: # check compartibility
                    if self.header[i][field]==header[field]:
                        raise ValueError('Cannot merge. Field: ' + field + ' is equal for multiple files.')

            # append header to current file
            self.header.append(header)

    def __add__(self, other):
        """add other Halos instance data to current instance"""
        self._add_header(other.header[0])
        self.files += other.files
        self.halos += other.halos
        self.h += other.h
        self.warnings = {}

    def __radd__(self, other):
        if other==0:
            return self
        else:
            return self.__add__(other)

    def __repr__(self):
        return '< ' + str(len(self.halos)) + ' halo collection >'

    def __eq__(self, other):
        return set([h.id for h in self.h]) == set([h.id for h in other.h])


class Halo(object):
    coord_type = COORDS

    def __init__(self, id, pos, particles):
        self.id = id
        self.pos = np.array(pos, dtype=Halo.coord_type).view(np.recarray)
        self.particles = particles  # record array of particles (see halos/helpers/helpers.py -> create_halo())
        self.particlesn = particles     # record array of particles transformed to new basis
        self.cov = np.identity(3, dtype=np.float32)   # covariance matrix of entire halo
        self.eig = np.identity(3, dtype=np.float32)    # eigenvectors of entire halo
        self.evals = np.empty(3, dtype=np.float32)       # eigenvalues of entire halo
        self.inner_cov = np.identity(3, dtype=np.float32)   # covariance matrix of half mass halo
        self.inner_eig = np.identity(3, dtype=np.float32)    # eigenvectors of half mass halo
        self.inner_evals = np.empty(3, dtype=np.float32)       # eigenvalues of half mass halo
        self.radii = np.empty(len(self.particles), dtype=np.float32)    # ellipsoidal radii of particles
        self.half_mass_radius = np.float32(0)                                  # ellipsoidal radius of hald mass sub-halo
        self.inner_R = np.array((0,0,0), dtype=Halo.coord_type)       # principal axis dimensions for inner halo
        self.fig = None
        self._center_flag = False                            # relative coords?

    def __repr__(self):
        return '< Halo: ' + str(self.id) + '; Size: ' + str(len(self.particles)) + ' particles >'

    def center_halo(self):
        """
        Iterating over each halo and subtracting its centre position
        from particle coordinates i.e. convert to local coordinates.
        """
        if not self._center_flag:
            self._center_flag = True
            self.particles.x -= self.pos.x
            self.particles.y -= self.pos.y
            self.particles.z -= self.pos.z

    def uncenter_halo(self):
        """
        Restore halo particles to absolute coordinates.
        """
        if self._center_flag:
            self._center_flag = False
            self.particles.x += self.pos.x
            self.particles.y += self.pos.y
            self.particles.z += self.pos.z

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
        self.inner_evals = self.evals
        self.eig = np.vstack((eigvecs[:, order[0]], eigvecs[:, order[1]], eigvecs[:, order[2]])).T
        #self.inner_eig = self.eig

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

        indices, _ = self.cut()
        self.inner_R = self.cleave(indices)       # calculate dimensions of inner halo

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
        return self.half_mass_radius

    def higher_order_fit(self, order=2):
        """
        Create sub-halo from particles inside half-mass radius, recompute covariance matrix, recalculate particles
        inside half-mass radius, repeat.
        :order number of times -1 to repeat
        """
        if order < 2:
            return
        for i in range(1, order):
            indices, _ = self.cut()
            # halo of only inner particles
            coords =  zip(self.particlesn.x[indices], self.particlesn.y[indices], self.particlesn.z[indices])
            h = Halo(id='subhalo', pos=(0,0,0), particles=np.array(coords, dtype=self.coord_type).view(np.recarray))
            # copy of self
            othercoords =  zip(self.particlesn.x, self.particlesn.y, self.particlesn.z)
            H = Halo(id='otherhalo', pos=(0,0,0), particles=np.array(othercoords, dtype=self.coord_type).view(np.recarray))

            h.get_covariance_matrix()              # calculate descriptors of inner halo
            h.get_eigenvectors()

            H.convert_basis(basis=h.eig)           # transform copy of self into new descriptors
            H.get_radii(evals=h.evals)
            indices, _ = H.cut()                   # find new 'inner' particles

            self.inner_cov = h.cov
            self.inner_evals = h.evals
            self.inner_eig = h.eig
            self.inner_R = H.inner_R
            self.radii[indices] = H.radii[indices]  # reassign subhalo radii to current halo's inner radii

        self.get_half_mass_radius()
        self.inner_R = h.cleave()

    def encapsulation(self, mode='cleave', transform=True):
        """
        Get a percentage of points that are not encapsulated by the ellipsoid, but should be.
        """
        if mode=='cleave':
            indices, _ = self.cut()
            t_matrix = np.identity(3) if not transform else np.linalg.inv(self.inner_eig)

            coords = np.array(zip(self.particlesn.x[indices], self.particlesn.y[indices], self.particlesn.z[indices])).T
            coords = np.dot(t_matrix, coords)

            ell_R = np.array([1./self.inner_R.x**2,  1./self.inner_R.y**2, 1./self.inner_R.z**2])
            coords2 = coords * coords
            ell_coords = np.dot(ell_R, coords2)
            r2 = np.einsum('ij->j', coords)
            total = np.sum(np.where(r2<=1, 1, 0))

            return format(total / float(len(indices)), '0.2f')

    def cut(self, fraction=1.0):
        """
        Returns indices of the first and second halves of particle/radius arrays corresponding to particles
        in and outside the half mass radius.
        :fraction fraction of particles to select from each cut
        """
        sortedindices = np.array(np.argsort(self.radii))
        l = int(len(sortedindices)/2)
        firsthalf = np.random.choice(sortedindices[:l], replace=False, size=int(l*fraction))
        secondhalf = np.random.choice(sortedindices[l:], replace=False, size=len(sortedindices)-int(l*fraction))
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

    def visualize(self, ellipsoids=False, mode='cleave', transform=True, fraction=1.0):
        """
        3D plot of particles. Particles within half mass radius are in red. Others are in blue.
        :fraction: fraction of particles to show (1.0 means 100%)
        :ellipsoids: whether to draw fitted surfaces
        :mode: whether to draw ellipsoids by largest dimension (cleave) or eigenvalues(eval)
        :transform: whether or not to transform inner ellipsoid to inner particle distribution
        """
        self.fig = plt.figure(self.id).add_subplot(111, projection='3d')
        firsthalf, secondhalf = self.cut(fraction)
        if self.half_mass_radius==0:
            colors = ('r','b')
        else:
            colors = ('b','b')
        self.fig.scatter(self.particlesn.x[firsthalf], self.particlesn.y[firsthalf], self.particlesn.z[firsthalf], colors[0])
        self.fig.scatter(self.particlesn.x[secondhalf], self.particlesn.y[secondhalf], self.particlesn.z[secondhalf], colors[1])
        if ellipsoids:
            self._draw_ellipsoids((firsthalf, secondhalf), mode, transform)
        plt.show()
        plt.close()

    def _draw_ellipsoids(self, indices, mode, transform=True):
        ax = self.fig

        transforms = [np.identity(3), np.identity(3)]

        if mode == 'cleave':            # ellipsoid radii determined on farthest points
            radii = [self.inner_R, self.cleave()] if transform else [self.cleave(indices[0]), self.cleave()]
        elif mode == 'eval':            # ellipsoid radii determined on distribution eigenvalues
            radii = [np.array((tuple(self.inner_evals * self.half_mass_radius / np.amax(self.inner_evals))), dtype=Halo.coord_type).view(np.recarray), \
                        np.array((tuple(self.evals * np.amax(self.radii) / np.amax(self.evals))), dtype=Halo.coord_type).view(np.recarray)]

        transforms = [self.inner_eig, np.identity(3)] if transform else transforms

        colors = ('r', 'c')
        alpha = 0.3

        u = np.linspace(0, 2 * np.pi, 100)              # generate ellipsoid surface data parametrically
        v = np.linspace(0, np.pi, 100)

        for r, c, t in zip(radii, colors, transforms):
            x = r.x * np.outer(np.cos(u), np.sin(v))
            y = r.y * np.outer(np.sin(u), np.sin(v))
            z = r.z * np.outer(np.ones_like(u), np.cos(v))
            tcoords = np.dot(t, np.array(zip(x, y, z)).T)        # ellipsoid coordinates
            ax.plot_surface(tcoords[0, :], tcoords[1, :], tcoords[2, :], rstride=4, cstride=4, color=c, alpha=alpha)
            alpha /= 2.0
        max_r = np.amax([radii[1].x, radii[1].y, radii[1].z])
        ax.set_xlim3d(-max_r, max_r)   # set axes limits
        ax.set_ylim3d(-max_r, max_r)
        ax.set_zlim3d(-max_r, max_r)

    def report(self):
        """
        Prints some stats about the halo
        """
        R = self.cleave()
        try:
            angles = np.arccos(np.einsum('ij,ij->j', self.inner_eig, self.eig))/np.pi
        except RuntimeWarning:
            angles = np.array([0,0,0])
        angles = [self._fmt(a) + 'pi' for a in angles]
        print '=' * 60
        print
        print '    Halo ID:\t\t\t' + str(self.id)
        print '    Particles:\t\t\t' + str(len(self.particles))
        print '    Half Mass Radius:\t\t' + str(self.half_mass_radius)
        print
        print '    Approx. Halo Dimensions:\t' + self._fmt(float(R.x)) + ' x ' + self._fmt(float(R.y)) + ' x '  + self._fmt(float(R.z))
        print '    Halo eigenvalues:\t\t' + self._fmt(self.evals[0]) + ', ' + self._fmt(self.evals[1]) + ', ' + self._fmt(self.evals[2])
        print
        print '    Approx inner dimensions:\t' + self._fmt(float(self.inner_R.x)) + ' x ' + self._fmt(float(self.inner_R.y)) + ' x '  + self._fmt(float(self.inner_R.z))
        print '    Inner eigenvalues:\t\t' + self._fmt(self.inner_evals[0]) + ', ' + self._fmt(self.inner_evals[1]) + ', ' + self._fmt(self.inner_evals[2])
        print '    Inner basis rotation:\t' + 'x->x ' + angles[0] + '; y->y ' + angles[1]  + '; z->z ' + angles[2]
        print
        print '=' * 60

    def _fmt(self, num, precision=2):
        """returns formated number according to precision
        """
        f_string = '{:.' + str(precision) + 'e}'
        return f_string.format(num)

    def __eq__(self, other):
        """checks if two Halo instances have same attributes
        :param other: instance of Halo object
        """
        return self.__dict__ == other.__dict__

    def intersect(self, other):
        """returns set of ids of particles in both Halo instances. Particles are
        identified by their IDs, not coordinates.
        :param other: instance of Halo object
        """
        a = set(self.particle.id)
        b = set(other.particle.id)
        return a & b

    def subtract(self, other):
        """removes particles contained in other from self. Particles are
        identified by their IDs, not coordinates.
        :param other: instance of Halo object
        """
        self.particles = self.particles[np.array([id in other.particles.id for \
                                        id in self.particles.id], dtype=bool)]

    def add(self, other):
        """add new particles to current Halo. Particles are identified by their
        IDs, not coordinates. Other particles might be in their local coordinate
        system which needs to be accounted for by uncenter_halo().
        :param other: a Halo instance
        """
        newparts = other.particles[np.array([id not in self.particles.id for \
                                    id in other.particles.id], dtype=bool)]
        np.resize(len(self.particles)+len(newparts))
        self.particles[len(self.particles):] = newparts

    def find_center(self):
        """take an arithmetic mean of particle coordinates to find the center
        position of halo
        """
        self.pos.x = np.mean(self.particles.x)
        self.pos.y = np.mean(self.particles.y)
        self.pos.z = np.mean(self.particles.z)
