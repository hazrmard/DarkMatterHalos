#!/usr/bin/env python

import sys
import bgc2
import numpy as np

mass_col = 0
pos_cols = (1,2,3)
vel_cols = (4,5,6)
halo_id_col = 0

dist_scale = 1.0e3
common_mass = 5.33423e5

good_halos = np.array([158801, 163119, 163651, 154976, 143362, 156963, 148512, 63791, 63985, 76713])
bad_halos = np.array([64077, 132222, 143444, 143452, 148568, 148569, 148571, 148572, 143446, 143454])


def main():
  for input_file in sys.argv[1:]:
    header, halos, particles = bgc2.read_bgc2(input_file)
    halos = np.asarray(halos)

    for index in range(len(halos)):
      halo_id = halos[index][halo_id_col]
      if (len(particles[index]) > 100) and ((halo_id == good_halos).any() or (halo_id == bad_halos).any()):
        #r_vir = halos[index][4] * dist_scale
        #halo_mass = halos[index][5]
        halo_pos = np.array([halos[index][6] * dist_scale, halos[index][7] * dist_scale, halos[index][8] * dist_scale])

        halo_particles = np.asarray(particles[index])
        pos = halo_particles[:,pos_cols[0]:pos_cols[0]+3] * dist_scale
        pos = pos - halo_pos
        mass = np.ones(len(halo_particles)) * common_mass

        print 'Using %d particles in halo %d.' % (halo_particles.shape[0], halo_id)

        if (halo_id == good_halos).any():
          outfile = 'good_halo_%d.dat' % (halo_id)
        elif (halo_id == bad_halos).any():
          outfile = 'bad_halo_%d.dat' % (halo_id)

        np.savetxt(outfile, np.column_stack((mass, pos)))



if __name__ == '__main__':
  main()

