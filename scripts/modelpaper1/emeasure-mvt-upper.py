import math
import scipy.ndimage as ndimage
import matplotlib.ticker as ticker
import numpy as np
from astropy.io import fits

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")
load_src("fmt", "../torchpack/formatting.py")
import fmt
load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData()
import lupa
lua = lupa.LuaRuntime(unpack_returned_tuples=True)

import torch
import hgspy

DPI = 300
figformat = 'png'
plot_size = 8
fontsize = 7
torch.set_font_sizes(fontsize)
cmap = hgspy.get_grey_cmap(isReversed=True)

snapshots = [10, 20, 30, 40, 50]

def addImage(col, row, grid, xrange):
	filename = mp1_data.getRadioDirname(row, 2, snapshots[col], 45, 5) + "/emeasure_ff.fits"

	ax = grid.grid[col][row]
	cbar_ax = grid.cgrid[col][row]

	hdu_list = fits.open(filename)
	image_data = hdu_list[0].data

	sig = 5
	FWHM = sig*hdu_list[0].header['PIXAS']*2.335
	image_data = ndimage.gaussian_filter(image_data, sigma=(sig, sig), order=0)

	nx = hdu_list[0].header['NAXIS1']
	ny = hdu_list[0].header['NAXIS2']

	dx = -hdu_list[0].header['CDELT1']*60.0*60.0
	dy = hdu_list[0].header['CDELT2']*60.0*60.0

	xmax = 0.5 * nx * dx
	xmin = -0.5 * nx * dx
	ymax = 0.5 * ny * dy
	ymin = -0.5 * ny * dy

	x = np.arange(xmin + (dx/2.0), xmax, dx)
	y = np.arange(ymin + (dy/2.0), ymax, dy)
	X, Y = np.meshgrid(x, y)
	max = hdu_list[0].header['MPIX']

	nlevels = 7
	levels = []

	for ilev in range(nlevels):
		levels.append(max/(math.sqrt(2.0)**(nlevels - ilev)))
	im = ax.imshow(image_data, extent=[xmin,xmax,ymin,ymax], origin='lower', cmap='gray_r')

	#cb = cbar_ax.colorbar(im, format=ticker.FuncFormatter(fmt))
	formatter = ticker.ScalarFormatter(useOffset=True, useMathText=True)
	formatter.set_powerlimits((0, 1))

	cb = grid.fig.colorbar(im, cax=cbar_ax, format=ticker.FuncFormatter(fmt.fmt), orientation='horizontal')
	#cb.ax.xaxis.set_ticks_position('top')
	xmax = cb.ax.get_xlim()[1]
	cb.ax.xaxis.set_ticks(np.arange(0, xmax, xmax/4.0))
	labels = cb.ax.get_xticklabels()
	labels[0] = ""
	labels[2] = ""
	cb.ax.set_xticklabels(labels)

	ax.contour(X, Y, image_data, levels, colors='black', linewidths=0.5)

	param_filename = mp1_data.getParamFilename(row, 2, snapshots[col])
	filestring = open(param_filename, 'r').read()

	table = lua.eval("{" + filestring + "}")

	physical_dx = table["Parameters"]["Grid"]["side_length"] / float(table["Parameters"]["Grid"]["no_cells_x"])
	side_length_y = table["Parameters"]["Grid"]["no_cells_y"] * physical_dx
	star_vert_phy = table["Parameters"]["Star"]["cell_position_y"] * physical_dx - 0.5 * side_length_y
	star_vert_ang = math.degrees(star_vert_phy / (1.5 * 1000 * 3.09e18)) * 60 * 60
	star_y = star_vert_ang / math.sqrt(2)
	cloud_vert_ang = star_vert_ang + math.degrees(0.35 / 1500.0) * 60 * 60
	cloud_y = cloud_vert_ang / math.sqrt(2)

	ax.scatter([0], [star_y], c='w', marker='*', linewidths=[0.4], zorder=2)
	ax.scatter([0], [cloud_y], c='w', marker='D', linewidths=[0.4], zorder=2)

	ax.set_xlim([-xrange[col][row]/2.0, xrange[col][row]/2.0])
	ax.set_ylim([-xrange[col][row]/2.0, xrange[col][row]/2.0])
	ax.xaxis.set_ticks(np.arange(-xrange[col][row]/2.0, (xrange[col][row]/2.0) + 0.0001, xrange[col][row]/4.0))
	ax.yaxis.set_ticks(np.arange(-(xrange[col][row]/2.0), (xrange[col][row]/2.0) + 0.0001, xrange[col][row]/4.0))

	ax.xaxis.get_major_ticks()[0].label1.set_visible(False)
	ax.xaxis.get_major_ticks()[2].label1.set_visible(False)
	ax.xaxis.get_major_ticks()[-1].label1.set_visible(False)

	ax.yaxis.get_major_ticks()[0].label1.set_visible(False)
	ax.yaxis.get_major_ticks()[2].label1.set_visible(False)
	ax.yaxis.get_major_ticks()[-1].label1.set_visible(False)

	if False:
		# Pixel angular size in upper left corner.
		ax.text(0.02, 0.96, '{:.2f}'.format(hdu_list[0].header['PIXAS']) + "''",
				fontsize=6, horizontalalignment='left', verticalalignment='top',
				rotation='horizontal', transform=ax.transAxes)

		# FWHM of guassian blur kernal in lower left corner.
		ax.text(0.02, 0.04, "FWHM=" + '{:.1f}'.format(FWHM) + "''",
				fontsize=6, horizontalalignment='left', verticalalignment='bottom',
				rotation='horizontal', transform=ax.transAxes)

		# FWHM diameter line to show scale
		ax.plot([(0.90)*xrange[col][row]/2.0 - FWHM, (0.90)*xrange[col][row]/2.0], [-(0.90)*xrange[col][row]/2.0, -(0.90)*xrange[col][row]/2.0])

		# Peak brightness in mJy/beam in upper right
		rbeam_rad = 0.5 * hdu_list[0].header['BMAJ'] * math.pi / 180.0
		pix_rad = abs(hdu_list[0].header['CDELT1']) * math.pi / 180.0

		intensity_filename = mp1_data.getRadioDirname(row, 2, snapshots[col], 45, 5) + "/intensity_pixel_ff.fits"
		hdu_list2 = fits.open(intensity_filename)

		peak = hdu_list2[0].header['MPIX'] * math.pi * rbeam_rad * rbeam_rad / (pix_rad * pix_rad)

		ax.text(0.98, 0.96, '{:.2f}'.format(peak) + "mJy/b",
				fontsize=6, horizontalalignment='right', verticalalignment='top',
				rotation='horizontal', transform=ax.transAxes)

	# Add column and row labels
	def fmt_mass(x):
		a = '{:.0f}'.format(x)
		return '{}'.format(a)

	if col == 0:
		ax.text(-0.22, 0.5, fmt.latexify("M = " + fmt.fmt_mass(mp1_data.masses[row]) + "\ \\mathrm{M_\\odot}"),
				fontsize=10, horizontalalignment='right', verticalalignment='center',
				rotation='vertical', transform=ax.transAxes)
	if row == 0:
		ax.text(0.5, 1.24, fmt.latexify("t = " + fmt.fmt_power(int(mp1_data.times[col] / 1000.0), '{:3.0f}', 0) + "\ \\mathrm{kyr}"),
				fontsize=10, horizontalalignment='center', verticalalignment='bottom',
				rotation='horizontal', transform=ax.transAxes)


hrats = [1, 1, 1, 1, 1]
vrats = [1, 1, 1, 1, 1]
fancy_grid = torch.FancyAxesGrid(hrats, vrats, border=(0.05, 0.01, 0.04, 0.02),
								 fig_width=plot_size, dpi=DPI, fig_format=figformat,
								 csize=0.010, cspace=0.004, cpad=0.03,
								 hspace=0.026, vspace=0.038)

xrange = []
xrange.append([0.040, 0.08, 0.18, 0.28, 0.40])
xrange.append([0.035, 0.10, 0.22, 0.34, 0.52])
xrange.append([0.038, 0.10, 0.24, 0.40, 0.65])
xrange.append([0.040, 0.10, 0.25, 0.48, 0.74])
xrange.append([0.040, 0.10, 0.25, 0.54, 0.80])

for i in range(5):
	for j in range(5):
		xrange[i][j] = xrange[i][j]*160.0

for i in range(5):
	for j in range(5):
		addImage(i, j, fancy_grid, xrange)

fancy_grid.update_cbar()

fancy_grid.save_plot("emeasure-mvt-upper.png")