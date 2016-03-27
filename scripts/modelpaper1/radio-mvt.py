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

import torch
import hgspy

DPI = 300
figformat = 'png'
plot_size = 10
torch.set_font_sizes(fontsize=8)
cmap = hgspy.get_grey_cmap(isReversed=True)

fileprefix = "data/model_paper1/set2/"

snapshots = [10, 20, 30, 40, 50]
param_index = range(19, 28, 1)

offgrid = np.genfromtxt(fileprefix + "offgrid_times", dtype=int)

masses = [6, 9, 12, 15, 20, 30, 40, 70, 120]
times = [2.0e4, 4.0e4, 6.0e4, 8e4, 10e4]

def latexify(str):
	return r'${}$'.format(str)

def fmt_nolatex(x, pos):
	if x >= -100 and x <= 100:
		a = '{:.2f}'.format(x)
		return '{}'.format(a)
	else:
		a, b = '{:.2e}'.format(x).split('e')
		b = int(b)
		return '{} \\times 10^{{{}}}'.format(a, b)

def fmt(x, pos):
	return latexify(fmt_nolatex(x, pos))

def getFilename(i, j):
	padi = "%03d" % (snapshots[i],)
	padj = "%02d" % (param_index[j],)
	return fileprefix + "data_" + padj + "/radio_" + padi + "/intensity_beam_ff.fits"

def addImage(col, row, grid, xrange):
	filename = getFilename(col, row)

	ax = grid.grid[col][row]
	cbar_ax = grid.cgrid[col][row]

	if offgrid[param_index[row] - 1] > snapshots[col]:
		hdu_list = fits.open(filename)
		image_data = hdu_list[0].data

		sig = 1 * math.pow(hdu_list[0].header['MPIX'] / 0.40, 0.14)
		FWHM = sig*hdu_list[0].header['PIXAS']*2.335
		image_data = ndimage.gaussian_filter(image_data, sigma=sig, order=0)

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

		im = ax.imshow(image_data, extent=[xmin,xmax,ymin,ymax], origin='lower', cmap='gray_r')

		#cb = cbar_ax.colorbar(im, format=ticker.FuncFormatter(fmt))
		formatter = ticker.ScalarFormatter(useOffset=True, useMathText=True)
		formatter.set_powerlimits((0, 1))

		cb = grid.fig.colorbar(im, cax=cbar_ax, format=ticker.FuncFormatter(fmt), orientation='horizontal')
		#cb.ax.xaxis.set_ticks_position('top')
		bmax = cb.ax.get_xlim()[1]
		cb.ax.xaxis.set_ticks(np.arange(0, bmax, bmax/4.0))
		labels = cb.ax.get_xticklabels()
		labels[0] = ""
		labels[2] = ""
		cb.ax.set_xticklabels(labels)

		vmax = hdu_list[0].header['MPIX']
		nlevels = 7
		levels = []

		for ilev in range(nlevels):
			levels.append(vmax/(math.sqrt(2.0)**(nlevels - ilev)))

		ax.contour(X, Y, image_data, levels, colors='black', linewidths=0.5)
	else:
		grid.set_visible(col, row, False)

		#cbar_ax.xaxis.set_ticks([])
		#cbar_ax.yaxis.set_ticks([])
		#cbar_ax.xaxis.set_visible(False)
		#cbar_ax.yaxis.set_visible(False)

		#ax.xaxis.set_ticks([])
		#ax.xaxis.set_ticks([])
		#ax.xaxis.set_visible(False)
		#ax.yaxis.set_visible(False)

	if grid.is_visible(col, row):
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

		peak = hdu_list[0].header['MPIX']

		ax.text(0.98, 0.96, '{:.2f}'.format(peak) + "mJy/b",
				fontsize=6, horizontalalignment='right', verticalalignment='top',
				rotation='horizontal', transform=ax.transAxes)

	# Add column and row labels
	def fmt_mass(x):
		a = '{:.0f}'.format(x)
		return '{}'.format(a)

	if col == 0:
		ax.text(-0.22, 0.5, latexify("M = " + fmt_mass(masses[row]) + "\ \\mathrm{M_\\odot}"),
				fontsize=10, horizontalalignment='right', verticalalignment='center',
				rotation='vertical', transform=ax.transAxes)
	if row == 0:
		ax.text(0.5, 1.24, latexify("t = " + fmt_nolatex(times[col], 0) + "\ \\mathrm{yrs}"),
				fontsize=10, horizontalalignment='center', verticalalignment='bottom',
				rotation='horizontal', transform=ax.transAxes)


hrats = [1, 1, 1, 1, 1]
vrats = [1, 1, 1, 1, 1, 1, 1, 1]
fancy_grid = torch.FancyAxesGrid(hrats, vrats, border=(0.05, 0.01, 0.03, 0.02),
								 fig_width=plot_size, dpi=DPI, fig_format=figformat,
								 csize=0.005, cspace=0.002, cpad=0.03,
								 hspace=0.026, vspace=0.025)

xrange = []
xrange.append([0.040, 0.08, 0.18, 0.28, 0.40, 0.52, 0.7, 1.0])
xrange.append([0.035, 0.10, 0.22, 0.34, 0.52, 0.80, 1.0, 1.0])
xrange.append([0.038, 0.10, 0.24, 0.40, 0.65, 0.90, 1.0, 1.0])
xrange.append([0.040, 0.10, 0.25, 0.48, 0.74, 1.00, 1.0, 1.0])
xrange.append([0.040, 0.10, 0.25, 0.54, 0.80, 1.00, 1.0, 1.0])

for i in range(5):
	for j in range(8):
		xrange[i][j] = xrange[i][j]*160.0

for i in range(5):
	for j in range(8):
		addImage(i, j, fancy_grid, xrange)

fancy_grid.update_cbar()

fancy_grid.save_plot("intensity-pix-mvt.png")