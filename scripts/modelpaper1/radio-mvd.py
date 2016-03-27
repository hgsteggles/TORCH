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

fileprefix = "data/model_paper1/set2/"

offgrid = np.fromfile(fileprefix + "offgrid_times", dtype=int, sep='\n')

masses = [6, 9, 12, 15, 20, 30, 40, 70, 120]
densities = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]

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


def getFilename(i):
	padi = "%02d" % ((i+1),)
	return fileprefix + "data_" + padi + "/radio_025/intensity_beam_ff.fits"

def addImage(index, grid, xrange):
	col = int(index / 9)
	row = index % 9

	filename = getFilename(index)

	ax = grid.grid[col][row]
	cbar_ax = grid.cgrid[col][row]


	if offgrid[index] > 25:
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

		cb = grid.fig.colorbar(im, cax=cbar_ax, format=ticker.FuncFormatter(fmt), orientation='horizontal')
		#cb.ax.xaxis.set_ticks_position('top')
		xmax = cb.ax.get_xlim()[1]
		cb.ax.xaxis.set_ticks(np.arange(0, xmax, xmax/4.0))
		labels = cb.ax.get_xticklabels()
		labels[0] = ""
		labels[2] = ""
		cb.ax.set_xticklabels(labels)

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
		ax.text(0.5, 1.24, latexify("n_\\star = " + fmt_nolatex(densities[col], 0) + "\ \\mathrm{cm^{-3}}"),
				fontsize=10, horizontalalignment='center', verticalalignment='bottom',
				rotation='horizontal', transform=ax.transAxes)


###	Data set up.
vs_types = []
vsminmax = []
color_maps = []

for i in range(45):
	vs_types.append(torch.VarType("emeasure", False))
	hdu_list = fits.open(getFilename(i))
	image_data = hdu_list[0].data
	vsminmax.append([0, hdu_list[0].header['MPIX']])
	color_maps.append(hgspy.get_grey_cmap())

plotparams = torch.PlotParams(None, vs_types, vsminmax, True,
							  'linear', (9,5), color_maps, tight=True, detail="all")

### Plotting
plotter = torch.Plotter(5, 9, plot_size, figformat, DPI)

#grid = plotter.getGrid(plotparams)
#plotter.modifyGrid(grid, True)

hrats = [1, 1, 1, 1, 1]
vrats = [1, 1, 1, 1, 1, 1, 1, 1]
fancy_grid = torch.FancyAxesGrid(hrats, vrats, border=(0.05, 0.01, 0.03, 0.02),
								 fig_width=plot_size, dpi=DPI, fig_format=figformat,
								 csize=0.005, cspace=0.002, cpad=0.03,
								 hspace=0.026, vspace=0.025)

xrange = []
xrange.append([15, 32, 62, 90, 134, 120, 140, 160])
xrange.append([10, 22, 48, 80, 110, 140, 140, 160])
xrange.append([6, 15, 34, 60, 90, 120, 150, 160])
xrange.append([4, 9, 24, 46, 66, 100, 130, 160])
xrange.append([2.52, 6, 16, 30, 54, 80, 110, 160])

for i in range(45):
	if (i + 1) % 9 != 0:
		addImage(i, fancy_grid, xrange)

fancy_grid.update_cbar()

fancy_grid.save_plot("intensity-pixel-mvd.png")