#!/usr/bin/env python

import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize as opt
import bgc2
from matplotlib.patches import Circle
from matplotlib import patheffects
from mpl_toolkits.axes_grid1 import ImageGrid
from scipy.ndimage.filters import gaussian_filter
from ipdb import set_trace



def main():
	'''
	ellipsoids.py
	
	This program reads in halo and particle data from bgc2 files and halo shape
	data from the Rockstar's ASCII output.  Halo shape information generated
	from Rockstar is used to fit halos with appropriately-rotated ellipsoidal
	shells, and the half-mass radius is found from the ellipsoidal radius of
	the (n/2)th particle.

	The halo IDs and resulting half-mass radii are saved to disk to be added as
	an additional column in the master halo database.  Optionally, plots are
	generated to demonstrate the ellipsoid halo fitting.

	Requires bgc2.py as a dependency.
	'''

	#opts, args = get_args(sys.argv[1:])
	#output_file, bgc2_files, ascii_files = parse_args(opts, args)
	ascii_files = [sys.argv[1]]
	bgc2_files = [sys.argv[2]]

	for (ascii_file, bgc2_file) in zip(ascii_files, bgc2_files):
		#  read in halo ID, shape data, etc. from Rockstar output
		ascii_header, ascii_data = read_files(ascii_file, header_line=0, rec_array=True)
		#  read in bgc2 files and make arrays of halo and particle data
		bgc2_header, halos, particles = bgc2.read_bgc2_numpy(bgc2_file)

		#  find array to sort halos by number of particles to work from the biggest down
		halo_indices = np.argsort(halos.npart)
		halo_indices = halo_indices[::-1]
		if start_halo > 0:
			halo_indices = halo_indices[start_halo:]

		#  loop through halos to work on one halo and particle list at a time
		for iteration, halo_index in enumerate(halo_indices):
			#  exit loop if max iteration reached
			if (max_iteration > 0) and (iteration >= max_iteration):
				break

			#  get data for current halo
			halo = np.array(halos[halo_index]).view(np.recarray)
			halo_particles = particles[halo_index]
			ascii_halo = ascii_data[ascii_data.id == halo.id]

			#  check for id match and duplicate halos
			if len(ascii_halo) != 1:
				print "Error:  found %d ASCII halo matches for halo ID %d." % (len(ascii_halo), ascii_halo.id[0])
				continue

			#  skip halos with fewer than specified number of particles
			if (npart_threshold > 0) and (halo.npart < npart_threshold):
				print "Skipping remaining halos with fewer than %d particles." % npart_threshold
				break

			#  convert Mpc to kpc for halo and particle positions
			print "Converting units to kpc..."
			halo.radius = halo.radius * dist_scale
			for pos in halo.x, halo.y, halo.z, halo_particles.x, halo_particles.y, halo_particles.z:
				pos[...] = pos * dist_scale

			#  make particle positions relative to halo center
			print "Making particle positions relative to halo center..."
			for particle_pos, halo_pos in zip([halo_particles.x, halo_particles.y, halo_particles.z], [halo.x, halo.y, halo.z]):
				particle_pos[...] = particle_pos - halo_pos

			#  optionally, generate simplified fake halo data for testing/debugging
			if test_fake_halo:
				ascii_halo, halo, halo_particles = make_fake_halo(ascii_halo, halo, halo_particles)

			#  convert particle cartesian coordinates to (spherical or ellipsoidal) radii
			print "Converting particle positions to spherical radii..."
			r_sphere = np.sqrt((halo_particles.x)**2 + (halo_particles.y)**2 + (halo_particles.z)**2)
			if method == 'sphere':
				r = r_sphere
			elif method == 'ellipsoid':
				print "Rotating eigenvalue matrix of axis ratios..."
				ratios = get_rotated_ratios_matrix(ascii_halo)
				#ratios = get_best_minor_rotation(ascii_halo.rvir, ratios, \
				ratios = get_best_minor_rotation(ascii_halo.Rs, ratios, \
				                                 np.array([ascii_halo.Ax[0], ascii_halo.Ay[0], ascii_halo.Az[0]]), \
				                                 np.column_stack((halo_particles.x, halo_particles.y, halo_particles.z)))
				print "Converting particle positions to ellipsoidal radii..."
				r = get_ellipsoid_r(ratios, np.column_stack((halo_particles.x, halo_particles.y, halo_particles.z)))

			#  find half-mass radius
			r_half_mass_sphere, r_half_mass_ell = get_half_mass_r(r, r_sphere)

			#  save result to array for later output to file
			#  !! todo -- add this

			#  make plots
			if (iteration < num_halos_to_plot) and (generate_testing_plots or generate_paper_plots):
				print "Making plot set %d for halo ID %d..." % (iteration, halo.id)
				make_plots(iteration, halo, ascii_halo, halo_particles, ratios, r, r_half_mass_ell)

		#  save results to file
		#  !! todo -- add this

	print 'Finished.'
	return



def read_files(files, header_line=None, comment_char='#', rec_array=False):
	header = None
	data = None
	if type(files) == str:
		files = [files]

	if header_line != None:
		with open(files[0], 'r') as fd:
			for line in range(header_line):
				fd.readline()
			header = fd.readline()
		if header[0] != comment_char:
			print "Header must start with a '%s'." % comment_char
			sys.exit(4)
		header = header[1:]
		header = header.split()

	for file in files:
		print "Reading file%s..." % (file)
		if data == None:
			if rec_array:
				data = np.genfromtxt(file, dtype=None, comments=comment_char, names=header, deletechars='[]/|')
				data = data.view(np.recarray)
			else:
				data = np.genfromtxt(file, dtype=None, comments=comment_char)
		else:
			if rec_array:
				data = np.append(data, np.genfromtxt(file, dtype=None, comments=comment_char, names=header, deletechars='[]/|'), axis=0)
				data = data.view(np.recarray)
			else:
				data = np.append(data, np.genfromtxt(file, dtype=None, comments=comment_char), axis=0)

	if header_line == None:
		return data
	else:
		return header, data



def get_rotated_ratios_matrix(ascii_halo):
	#  get rotation angles from A vector (and make them numbers instead of np arrays)
	theta_z = atan(ascii_halo.Ay, ascii_halo.Ax)[0]
	theta_y = atan(-ascii_halo.Az, ascii_halo.Ax)[0]
	#theta_x = atan(ascii_halo.Az, ascii_halo.Ay)[0]

	#  form a diagonal matrix of the inverse-squared axis ratios
	ratios = np.diag([1.0, 1.0 / (ascii_halo.b_to_a)**2, 1.0 / (ascii_halo.c_to_a)**2])

	#  rotate axis ratio matrix about z-axis  -->  X_rot = R^T * X * R
	ratios = z_rotation_matrix(theta_z).dot(ratios).dot(z_rotation_matrix(theta_z).T)
	#ratios = z_rotation_matrix(theta_z).T.dot(ratios).dot(z_rotation_matrix(theta_z))

	#  rotate axis ratio matrix about y-axis  -->  X_rot = R^T * X * R
	ratios = y_rotation_matrix(theta_y).dot(ratios).dot(y_rotation_matrix(theta_y).T)
	#ratios = y_rotation_matrix(theta_y).T.dot(ratios).dot(y_rotation_matrix(theta_y))

	#  rotate axis ratio matrix about x-axis  -->  X_rot = R^T * X * R
	#ratios = x_rotation_matrix(theta_x).dot(ratios).dot(x_rotation_matrix(theta_x).T)
	#ratios = x_rotation_matrix(theta_x).T.dot(ratios).dot(x_rotation_matrix(theta_x))

	return ratios



def get_best_minor_rotation(radius, ratios, A, pos):
	theta = get_best_theta(radius, ratios, A, pos)
	ratios = rotate_ratios_about_axis(theta, A, ratios)
	return ratios



def get_best_theta(radius, ratios, A, pos):
	print 'Npart = ', len(pos)
	theta, fval, ierr, iters = opt.fminbound(count_npart_outside_radius, \
	                                         0.0, np.pi, \
	                                         args=([A, ratios, radius, pos]), \
	                                         xtol=1.e-2, maxfun=100, full_output=True, disp=3)
	if ierr == 0:
		return theta
	elif ierr == 1:
		print "Error in finding best minor rotation:  maximum number of function calls reached."
		return np.nan
	else:
		print "Unknown error flag received from fmim_bound().  Exiting..."
		sys.exit()



def count_npart_outside_radius(theta, A, ratios, radius, pos):
	ratios = rotate_ratios_about_axis(theta, A, ratios)
	r_ell = get_ellipsoid_r(ratios, pos)
	mask = (r_ell > radius)
	npart = mask.sum()
	return npart



def rotate_ratios_about_axis(theta, A, ratios):
	rotation = axis_angle_rotation_matrix(A[0], A[1], A[2], theta)
	ratios = rotation.dot(ratios).dot(rotation.T)
	return ratios



def get_ellipsoid_r(ratios, pos):
	#  convert particle cartesian coordinates to "ellipsoidal" radii
	#  using einstein notation here (sorry), since it's much faster and np.dot/np.tensordot do not behave as expected
	#  we want the element-wise dot product  ->  `pos.dot(ratios).dot(pos.T)`, but for each particle individually
	#  `temp = ratios.dot(pos.T)` works as expected, but `pos.dot(temp)` yields a (N, N) array, but we want (N, )
	#  np.einsum follows einstein notation of summation over repeated indices

	#  temporarily store results for `ratios <dot> pos`
	#  should be equivalent to `ratios.dot(pos.T)`
	tempdot = np.einsum('ij, jk -> ik', ratios, pos.T, order='C')

	#  now find `pos <dot> tempdot` for each row in pos
	#  should be equivalent to `pos[:,0]*tempdot[0,:] + pos[:,1]*tempdot[1,:] + pos[:,2]*tempdot[2,:]`
	r_ell = np.einsum('ij, ji -> i', pos, tempdot, order='C')

	#  take square root and get rid of the extra array dimension to make 1-D
	r_ell = np.sqrt(r_ell)
	r_ell = np.squeeze(r_ell)

	return r_ell



def get_half_mass_r(r, r_real=None):
	#  set r_real to r if using spherical radii and not already done
	if r_real == None:
		r_real = r

	#  find (n/2)th particle(s) and corresponding half-mass radius
	sort_indices = np.argsort(r)
	if len(r) % 2 != 0:
		#  if odd number of particles, simply find the radius of the middle particle
		half_index = sort_indices[len(sort_indices) / 2]
		r_half_mass_sphere = r_real[half_index]
		r_half_mass_ell = r[half_index]
	elif len(r) % 2 == 0:
		#  if even number of particles, find two middle particles and radius halfway between them
		half_index1 = sort_indices[len(sort_indices) / 2 - 1]
		half_index2 = sort_indices[len(sort_indices) / 2]
		r_half_mass_sphere = (r_real[half_index1] + r_real[half_index2]) / 2.
		r_half_mass_ell = (r[half_index1] + r[half_index2]) / 2.

	return r_half_mass_sphere, r_half_mass_ell



def x_rotation_matrix(theta):
	return np.array([[1., 0.,           0.], \
	                 [0., cos(theta), -(sin(theta))], \
	                 [0., sin(theta),   cos(theta)]])



def y_rotation_matrix(theta):
	return np.array([[  cos(theta),  0., sin(theta)], \
	                 [  0.,          1., 0.], \
	                 [-(sin(theta)), 0., cos(theta)]])



def z_rotation_matrix(theta):
	return np.array([[cos(theta), -(sin(theta)), 0.], \
	                 [sin(theta),   cos(theta),  0.], \
	                 [0.,           0.,          1.]])



def axis_angle_rotation_matrix(x, y, z, theta):
	if np.sqrt(x**2 + y**2 + z**2) != 1.:
		norm = np.sqrt(x**2 + y**2 + z**2)
		x, y, z = x / norm, y / norm, z / norm
	return np.array([[cos(theta) + x * x * (1. - cos(theta)),     x * y * (1. - cos(theta)) - z * sin(theta), x * z * (1. - cos(theta)) + y * sin(theta)], \
	                 [y * x * (1. - cos(theta)) + z * sin(theta), cos(theta) + y * y * (1. - cos(theta)),     y * z * (1. - cos(theta)) - x * sin(theta)], \
	                 [z * x * (1. - cos(theta)) - y * sin(theta), z * y * (1. - cos(theta)) + x * sin(theta), cos(theta) + z * z * (1. - cos(theta))]])



def sin(theta):
	return np.sin(theta)



def cos(theta):
	return np.cos(theta)



def tan(theta):
	return np.tan(theta)



def atan(x1, x2):
	return np.arctan2(x1, x2)



def make_plots(iteration, halo, ascii_halo, particles, ratios, r, r_half_mass):
	id_string = "%.2d--id%d" % (iteration, halo.id)
	if generate_testing_plots:
		ellipse_ring_fig = make_ellipse_ring_plot(id_string, particles, r, r_half_mass)
		#projectinos_fig  = make_projections_plot(id_string, particles, r_half_mass, halo.radius)
		in_out_fig       = make_in_out_plot(id_string, halo, particles, r, r_half_mass)
	if generate_paper_plots:
		proj_in_out_grid = make_proj_in_out_grid_plot(id_string, halo, ascii_halo, particles, r, r_half_mass)

	return 0



def make_ellipse_ring_plot(id_string, particles, r, r_half_mass):
	fig = plt.figure(figsize = (9.0, 6.0))
	ax = fig.add_subplot(111, aspect='equal')

	ax = draw_ellipse_ring(ax, particles, r, r_half_mass)

	fig.tight_layout()
	plot_name = "%s%s%s%s" % (plot_base, 'ellipse_ring-', id_string, plot_ext)
	plt.savefig(plot_name, bbox_inches='tight')

	return fig



def draw_ellipse_ring(ax, particles, r, r_half_mass):
	mask = get_slice_mask(particles.z)
	particles = particles[mask]
	r = r[mask]

	mask = (np.abs(r - r_half_mass) / r_half_mass <= 0.05)
	particles = particles[mask]

	ax.plot(particles.x, particles.y, linestyle='', marker='.', markersize=2, markeredgecolor='blue')
	ax.add_patch(Circle((0., 0.), r_half_mass, fc="None", ec="black", lw=1))

	return ax



def get_slice_mask(z):
	mask = (np.abs(z) <= z.max() * slice_fraction)
	return mask




def make_projections_plot(id_string, particles, r_half_mass, r_vir):
	fig = plt.figure(figsize = (9.0, 6.0))
	ax = fig.add_subplot(111, aspect='equal')

	plot_lim = np.max((particles.x.max(), particles.y.max()))
	plot_lim = plot_lim + plot_lim * 0.1
	ax = draw_projection(ax, particles.x, particles.y, plot_lim)
	ax.add_patch(Circle((0., 0.), r_half_mass, fc="None", ec="black", lw=1))
	ax.add_patch(Circle((0., 0.), r_vir, fc="None", ec="black", lw=1))

	fig.tight_layout()
	plot_name = "%s%s%s%s" % (plot_base, 'projections-', id_string, plot_ext)
	plt.savefig(plot_name, bbox_inches='tight')

	return fig



def draw_projection(ax, x, y, plot_lim, hx = None, hy = None, r = None):
	limits = [[-plot_lim, plot_lim], [-plot_lim, plot_lim]]
	z, xedges, yedges = np.histogram2d(x, y, bins=npixels, range=limits)
	if log_scale_projections:
		z[z<1.0] = 0.5
		plot_norm = mpl.colors.LogNorm(vmin = 1, vmax = z.max(), clip=True)
	else:
		plot_norm = None
	if extra_smoothing:
		z = gaussian_filter(z, smoothing_radius)
	im = ax.imshow(z.T, extent=(-plot_lim, plot_lim, -plot_lim, plot_lim), \
			interpolation='gaussian', origin='lower', cmap=colormap, norm=plot_norm)
	ax.locator_params(nbins=6)
	if draw_circle and hx != None and hy != None and r != None:
		ax.add_patch(Circle((hx, hy), r, fc="None", ec="black", lw=1))
	if draw_contours:
		x_midpoints = (xedges[:-1] + xedges[1:]) / 2.0
		y_midpoints = (yedges[:-1] + yedges[1:]) / 2.0
		X, Y = np.meshgrid(x_midpoints, y_midpoints)
		ax.contour(X, Y, z.T, 2, colors='black', linewidths=4)
		ax.contour(X, Y, z.T, 2, colors='white', linewidths=2)
	if label_colorbar:
		if log_scale_projections:
			log_format = mpl.ticker.LogFormatterMathtext(10, labelOnlyBase=False)
			ax.cax.colorbar(im, format=log_format)
		else:
			ax.cax.colorbar(im)
	else:
		bar = ax.cax.colorbar(im, ticks=[])
		bar.ax.set_yticklabels([])
		#plt.setp(bar.ax.get_yticklabels(), visible=False)

	return ax



def make_in_out_plot(id_string, halo, particles, r, r_half_mass):
	fig = plt.figure(figsize = (9.0, 6.0))
	ax = fig.add_subplot(111, aspect='equal')

	plot_lim = halo.radius
	plot_lim = plot_lim + plot_lim * 0.1

	mask = get_slice_mask(particles.z)
	particles = particles[mask]
	r = r[mask]

	ax = draw_in_out_particles(ax, particles.x, particles.y, r, r_half_mass, plot_lim)
	ax.add_patch(Circle((0., 0.), halo.radius, fc="None", ec="black", lw=1))

	fig.tight_layout()
	plot_name = "%s%s%s%s" % (plot_base, 'in_out-', id_string, plot_ext)
	plt.savefig(plot_name, bbox_inches='tight')

	return fig



def draw_in_out_particles(ax, x, y, r, r_half_mass, plot_lim):
	mask = (r < r_half_mass)
	x_in = x[mask]
	y_in = y[mask]

	mask = (r > r_half_mass)
	x_out = x[mask]
	y_out = y[mask]

	ax.plot(x_in, y_in, linestyle='', marker='.', markersize=1, markeredgecolor='blue')
	ax.plot(x_out, y_out, linestyle='', marker='.', markersize=1, markeredgecolor='red')

	ax.set_xlim([-plot_lim, plot_lim])
	ax.set_ylim([-plot_lim, plot_lim])

	return ax



def make_proj_in_out_grid_plot(id_string, halo, ascii_halo, particles, r, r_half_mass):
	fig = plt.figure(figsize = (9.0, 6.0))

	plot_lim = halo.radius
	plot_lim = plot_lim + plot_lim * 0.1

	if label_projections:
		#ax = fig.add_subplot(211, aspect=2.0/3.2)
		ax = fig.add_subplot(111)
		ax = hide_axes(ax)
		ax.set_xlabel(proj_xlabel)
		ax.set_ylabel(proj_ylabel)

	fig = make_projections(fig, 211, ascii_halo, particles, halo.radius, plot_lim)
	fig = make_in_out_panels(fig, 212, particles, r, r_half_mass, halo.radius, plot_lim)

	fig.tight_layout()
	plot_name = "%s%s%s%s" % (plot_base, 'proj_grid-', id_string, plot_ext)
	plt.savefig(plot_name, bbox_inches='tight')

	return fig



def make_projections(fig, position, ascii_halo, particles, r_vir, plot_lim):
	grid = ImageGrid(fig, position, nrows_ncols=(1,3), axes_pad=0.12, cbar_mode='single')
	for i, (x, y, z, Ax, Ay, Az) in enumerate(zip( \
	        (particles.x, particles.x, particles.y), \
	        (particles.y, particles.z, particles.z), \
	        (particles.z, particles.y, particles.x), \
	        (ascii_halo.Ax, ascii_halo.Ax, ascii_halo.Ay), \
	        (ascii_halo.Ay, ascii_halo.Az, ascii_halo.Az), \
	        (ascii_halo.Az, ascii_halo.Ay, ascii_halo.Ax))):
		ax = grid[i]

		if slice_projections:
			mask = get_slice_mask(z)
			x = x[mask]
			y = y[mask]

		ax = draw_projection(ax, x, y, plot_lim)
		ax.add_patch(Circle((0., 0.), r_vir, fc="None", ec="black", lw=1))
		ax = draw_ellipsoid_vector(ax, Ax, Ay, Az, r_vir)

		if print_text:
			if i == 0:
				ax.text(0.95, 0.88, 'XY', color='black', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=2, foreground='white')])
			if i == 1:
				ax.text(0.95, 0.88, 'XZ', color='black', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=2, foreground='white')])
			if i == 2:
				ax.text(0.95, 0.88, 'YZ', color='black', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=2, foreground='white')])
	return fig



def make_in_out_panels(fig, position, particles, r, r_half_mass, r_vir, plot_lim):
	grid = ImageGrid(fig, position, nrows_ncols=(1,3), axes_pad=0.12, cbar_mode='single')
	for i, (x, y, z) in enumerate(zip( \
	        (particles.x, particles.x, particles.y), \
	        (particles.y, particles.z, particles.z), \
	        (particles.z, particles.y, particles.x))):
		ax = grid[i]

		mask = get_slice_mask(z)
		x_slice = x[mask]
		y_slice = y[mask]
		r_slice = r[mask]

		ax = draw_in_out_particles(ax, x_slice, y_slice, r_slice, r_half_mass, plot_lim)
		ax.add_patch(Circle((0., 0.), r_vir, fc="None", ec="black", lw=1))

		#bar = ax.cax.colorbar(im, ticks=[])
		bar = mpl.colorbar.ColorbarBase(ax.cax, cmap='seismic', ticks=[])
		bar.ax.set_yticklabels([])

		if print_text:
			if i == 0:
				ax.text(0.95, 0.88, 'XY', color='black', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=2, foreground='white')])
			if i == 1:
				ax.text(0.95, 0.88, 'XZ', color='black', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=2, foreground='white')])
			if i == 2:
				ax.text(0.95, 0.88, 'YZ', color='black', horizontalalignment='right', verticalalignment='center', transform=ax.transAxes, path_effects=[patheffects.withStroke(linewidth=2, foreground='white')])
	return fig



def draw_ellipsoid_vector(ax, Ax, Ay, Az, r_vir):
	scale = r_vir / np.sqrt(Ax**2 + Ay**2 + Az**2)
	Ax = Ax * scale
	Ay = Ay * scale
	if Ax != 0 or Ay != 0:
		ax.arrow(0., 0., Ax[0], Ay[0], length_includes_head=True, linewidth=1, width=r_vir/20., facecolor='0.75', edgecolor='black')
	return ax



def hide_axes(ax):
	ax.spines['top'].set_color('none')
	ax.spines['bottom'].set_color('none')
	ax.spines['left'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')
	return ax



def make_fake_halo(ascii_halo, halo, particles):
	r_vir = 50.
	halo.radius = r_vir
	ascii_halo.rvir = r_vir
	ascii_halo.b_to_a = 0.5
	ascii_halo.c_to_a = 0.25
	ascii_halo.Ax = r_vir
	ascii_halo.Ay = 0.0
	ascii_halo.Az = 0.0
	particles.x = r_vir * (np.random.random(len(particles)) * 2. - 1.)
	particles.y = r_vir * (np.random.random(len(particles)) * 1. - 0.5)
	particles.z = r_vir * (np.random.random(len(particles)) * 0.5 - 0.25)

	rotation_matrices = []
	for rotation in rotate_fake_order:
		if rotation == 'x':
			rotation_matrices.append(x_rotation_matrix(rotate_fake_x_angle))
		if rotation == 'y':
			rotation_matrices.append(y_rotation_matrix(rotate_fake_y_angle))
		if rotation == 'z':
			rotation_matrices.append(z_rotation_matrix(rotate_fake_z_angle))

	for rotation_matrix in rotation_matrices:
		particles.x, particles.y, particles.z, ascii_halo.Ax, ascii_halo.Ay, ascii_halo.Az = rotate_fake_halo( \
		             rotation_matrix, \
		             particles.x, particles.y, particles.z, \
		             ascii_halo.Ax, ascii_halo.Ay, ascii_halo.Az)

	return ascii_halo, halo, particles



def rotate_fake_halo(rotation_matrix, x, y, z, Ax, Ay, Az):
	pos = np.column_stack((x, y, z))
	pos = rotation_matrix.dot(pos.T).T
	A = np.array([Ax, Ay, Az])
	A = rotation_matrix.dot(A)
	return pos[:,0], pos[:,1], pos[:,2], A[0], A[1], A[2]




def add_white_to_colormap(orig_map, num):
	from matplotlib import cm
	temp_cmap = cm.get_cmap(orig_map, num)
	vals = temp_cmap(np.arange(num))
	nfade = num / 7
	vals[:nfade,0] = np.linspace(1., vals[nfade-1,0], nfade)
	vals[:nfade,1] = np.linspace(1., vals[nfade-1,1], nfade)
	vals[:nfade,2] = np.linspace(1., vals[nfade-1,2], nfade)
	newcmap = mpl.colors.LinearSegmentedColormap.from_list("custom_1", vals)
	return newcmap



colormap = add_white_to_colormap('rainbow', 30)
plot_dest_type = 'paper'
if plot_dest_type == 'paper':
	mpl.rcParams['font.family'] = 'serif'
	mpl.rcParams['font.size'] = 16
	mpl.rcParams['axes.linewidth'] = 3
	mpl.rcParams['lines.linewidth'] = 4
	mpl.rcParams['patch.linewidth'] = 4
	mpl.rcParams['xtick.major.width'] = 3
	mpl.rcParams['ytick.major.width'] = 3
	mpl.rcParams['xtick.major.size'] = 8
	mpl.rcParams['ytick.major.size'] = 8
	mpl.rcParams['xtick.minor.width'] = 2
	mpl.rcParams['ytick.minor.width'] = 2
	mpl.rcParams['xtick.minor.size'] = 4
	mpl.rcParams['ytick.minor.size'] = 4



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#	user-settable control parameters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
test_fake_halo = False			# generate an idealized fake halo for testing
if test_fake_halo:
	#rotate_fake_order = ['x', 'y', 'z']
	rotate_fake_order = ['x', 'z']
	rotate_fake_x_angle = np.pi / 6.
	rotate_fake_y_angle = -np.pi / 3.
	rotate_fake_z_angle = -np.pi / 3.
start_halo = 0					# first halo to analyze
max_iteration = 10				# number of halos to analyze
#max_iteration = None			# number of halos to analyze
num_halos_to_plot = 10			# max number of halos to make plots for
npart_threshold = 100			# minimum number of particles per halo
dist_scale = 1.e3				# convert Mpc to kpc
#method = 'sphere'				# use spherical shells for finding half-mass radius
method = 'ellipsoid'			# use ellipsoidal shells for finding half-mass radius
generate_testing_plots = True	# output plots for testing purposes
generate_paper_plots = True		# output plots for use in a paper
plot_base = 'plots/'			# string to prepend to plot file path
plot_ext = '.eps'				# string to append to plot file path
slice_fraction = 0.05
npixels = 250
smoothing_radius = 0.9
slice_projections = True
label_colorbar = False
draw_circle = True
draw_contours = False
log_scale_projections = True
extra_smoothing = True
print_text = True
label_projections = True
proj_xlabel = r'Position (kpc h$^{-1}$)'
proj_ylabel = r'Position (kpc h$^{-1}$)'
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__ == '__main__':
	main()
