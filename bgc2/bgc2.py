#!/usr/bin/env python

import sys


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
	for i in range(ngroups):
		group.append(list(struct.unpack(groupformat, group_bin[i*groupsize:(i+1)*groupsize])))

	# Particle stuff
	particlestart = headersize + groupoffset + ngroups*groupsize + particleoffset
	particle_bin = bin_string[particlestart:]
	particle = []
	p_start = 0
	for i in range(ngroups):
		npart = group[i][2]
		particle.append([])
		for j in range(npart):
			particle[i].append(list(struct.unpack(particleformat, particle_bin[p_start:p_start+particlesize])))
			p_start += particlesize
		p_start += particleoffset

	print "Finished parsing bgc2 file"
	return header, group, particle



def read_bgc2_numpy(filename):
	import numpy as np

	offset = 4
	groupoffset = 8
	particleoffset = 8

	headersize = 1024
	groupsize = 4*8 + 10*4
	particlesize = 1*8 + 6*4

	dt_header = np.dtype([('magic', np.uint64), \
	                      ('version', np.int64), \
	                      ('num_files', np.int64), \
	                      ('file_id', np.int64), \
	                      ('snapshot', np.int64), \
	                      ('format_group_data', np.int64), \
	                      ('format_part_data', np.int64), \
	                      ('group_type', np.int64), \
	                      ('ngroups', np.int64), \
	                      ('ngroups_total', np.int64), \
	                      ('npart', np.int64), \
	                      ('npart_total', np.int64), \
	                      ('npart_orig', np.int64), \
	                      ('max_npart', np.int64), \
	                      ('max_npart_total', np.int64), \
	                      ('min_group_part', np.int64), \
	                      ('valid_part_ids', np.int64), \
	                      ('linkinglength', np.float64), \
	                      ('overdensity', np.float64), \
	                      ('time', np.float64), \
	                      ('redshift', np.float64), \
	                      ('box_size', np.float64), \
	                      ('box_min_x', np.float64), \
	                      ('box_min_y', np.float64), \
	                      ('box_min_z', np.float64), \
	                      ('bounds_xmin', np.float64), \
	                      ('bounds_xmax', np.float64), \
	                      ('bounds_ymin', np.float64), \
	                      ('bounds_ymax', np.float64), \
	                      ('bounds_zmin', np.float64), \
	                      ('bounds_zmax', np.float64), \
	                      ('part_mass', np.float64), \
	                      ('Omega0', np.float64), \
	                      ('OmegaLambda', np.float64), \
	                      ('Hubble0', np.float64), \
	                      ('GravConst', np.float64)])

	dt_groups = np.dtype([('id', np.int64), \
	                      ('parent_id', np.int64), \
	                      ('npart', np.uint64), \
	                      ('npart_self', np.uint64), \
	                      ('radius', np.float32), \
	                      ('mass', np.float32), \
	                      ('x', np.float32), \
	                      ('y', np.float32), \
	                      ('z', np.float32), \
	                      ('vx', np.float32), \
	                      ('vy', np.float32), \
	                      ('vz', np.float32), \
	                      ('vmax', np.float32), \
	                      ('rvmax', np.float32)])

	dt_particles = np.dtype([('id', np.int64), \
	                         ('x', np.float32), \
	                         ('y', np.float32), \
	                         ('z', np.float32), \
	                         ('vx', np.float32), \
	                         ('vy', np.float32), \
	                         ('vz', np.float32)])

	print "Reading "+filename+"..."
	with open(filename, 'rb') as fd:
		fd.seek(offset, 0)

		# Header stuff
		header = np.rec.fromfile(fd, dtype=dt_header, shape=1)
		header = header[0]
		print 'This is bgc2 file %d of %d.' % (header.file_id, header.num_files)
		print 'Redshift = ', header.redshift
		print 'Number of halos = ', header.ngroups

		# Group/halo stuff
		fd.seek(offset + headersize + groupoffset, 0)
		groups = np.rec.fromfile(fd, dtype=dt_groups, shape=header.ngroups)

		# Particle stuff
		fd.seek(particleoffset, 1)
		particles = []
		for i in range(header.ngroups):
			particles.append(np.rec.fromfile(fd, dtype=dt_particles, shape=groups[i].npart))
			fd.seek(particleoffset, 1)

	print "Finished reading bgc2 file."
	return header, groups, particles



def main():
	header, group, particle = read_bgc2(sys.argv[1])

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
