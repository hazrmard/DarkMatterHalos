#!/usr/bin/env python

import sys
from .. import config
import numpy as np


def read_bgc2(filename):
	import struct

	offset = 4
	groupoffset = 8
	particleoffset = 8

	headersize = 1024
	groupsize = 4*8 + 10*4
	particlesize = 1*8 + 6*4

	headerformat = '=Q 16q 19d'
	groupformat = '=2q 2Q 10f'
	particleformat = '=q 6f'

	print "Reading "+filename+"..."
	fd = open(filename, 'rb')
	bin_string = fd.read()
	fd.close()
	print "Finished reading file."
	bin_string = bin_string[offset:]

	# Header stuff
	header_bin = bin_string[:headersize]
	header_pad = headersize - 36*8
	header = list(struct.unpack(headerformat, header_bin[:-header_pad]))

	# Group stuff
	ngroups = header[8]
	print 'ngroups = ', ngroups
	groupstart = headersize + groupoffset
	groupend = groupstart + ngroups*groupsize
	group_bin = bin_string[groupstart:groupend]
	group = []
	for i in xrange(ngroups):
		group.append(list(struct.unpack(groupformat, group_bin[i*groupsize:(i+1)*groupsize])))

	# Particle stuff
	particlestart = headersize + groupoffset + ngroups*groupsize + particleoffset
	particle_bin = bin_string[particlestart:]
	particle = []
	p_start = 0
	for i in xrange(ngroups):
		npart = group[i][2]
		particle.append([])
		for j in range(npart):
			particle[i].append(list(struct.unpack(particleformat, particle_bin[p_start:p_start+particlesize])))
			p_start += particlesize
		p_start += particleoffset

	print "Finished parsing bgc2 file"
	return header, group, particle



def read_bgc2_numpy(filename, level=2, sieve=None, onlyid=False):
	"""read halo data from bgc2 file.
	:filename - path to file
	:level - 0:only read header, 1:only read halo meta deta, 2: read particle data
	:sieve - set of halo ids to keep. If None, all ids are kept.
	:onlyid - only pass on ids for groups and particles if True. Usefule for preserving memory.
	"""
	import numpy as np

	offset = 4
	groupoffset = 8
	particleoffset = 8

	headersize = 1024
	groupsize = 4*8 + 10*4
	particlesize = 1*8 + 6*4

	header = None
	groups = None
	particles = None
	if sieve is not None:
		sieve = set(sieve)

	#print "Reading "+filename+"..."
	with open(filename, 'rb') as fd:
		fd.seek(offset, 0)

		if level>=0:
			# Header stuff
			header = np.rec.fromfile(fd, dtype=config.DT_HEADER, shape=1)
			header = header[0]


		if level>=1:
			# Group/halo stuff
			fd.seek(offset + headersize + groupoffset, 0)
			groups = np.rec.fromfile(fd, dtype=config.DT_GROUPS, shape=header.ngroups)
			if sieve is not None:	# set a filter flag list for groups
				temp_groups = np.array([g in sieve for g in groups.id], dtype=np.bool)


		if level>=2:
			# Particle stuff
			fd.seek(particleoffset, 1)
			particles = []
			if sieve is None:
				for i in range(header.ngroups):
					particles.append(np.rec.fromfile(fd, dtype=config.DT_PARTICLES, shape=groups[i].npart))
					fd.seek(particleoffset, 1)
			else:	# use the filter flag list to include or skip group particles
				for i in range(header.ngroups):
					# if not isinstance(temp_groups[i], (int,long)):
					if temp_groups[i]==1:	# needs to be included
						record = np.rec.fromfile(fd, dtype=config.DT_PARTICLES, shape=groups[i].npart)
						particles.append(record)
						fd.seek(particleoffset, 1)
					else:	# else num of particles in unneeded group so can be skipped
						fd.seek(particleoffset + groups[i].npart*particlesize, 1)

				groups = groups[temp_groups]	# finally filter groups as well using filter flag list


	if onlyid:
		particles = [np.array(p.id.copy(), dtype=np.dtype([('id',np.float64)])).view(np.recarray) for p in particles]
		groups = [g.id for g in groups]

	#print "Finished reading bgc2 file."
	return header, groups, particles



def main():
	header, group, particle = read_bgc2_numpy(sys.argv[1])

	print 'Header contents:'
	for value in header:
		print value
	print

	print 'Group[0] contents:'
	for value in group[0]:
		print value
	print

	print 'Particles in group[0]:'
	for part in particle[0]:
		print part
	print

	print 'Group[1] contents:'
	for value in group[1]:
		print value
	print

	print 'Particles in group[1]:'
	for part in particle[1]:
		print part





if __name__ == '__main__':
	main()
