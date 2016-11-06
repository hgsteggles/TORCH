import math
import scipy.ndimage as ndimage
import matplotlib.ticker as ticker
import numpy as np
from astropy.io import fits

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
import torch
load_src("hgspy", "../torchpack/hgspy.py")
import hgspy
load_src("fmt", "../torchpack/formatting.py")
import fmt
load_src("mp1", "../torchpack/modelpaper1.py")
import mp1

DPI = 300
figformat = 'png'
plot_size = 8
torch.set_font_sizes(fontsize=8)

angles = [0, 30, 45, 60, 90]

mp1_data = mp1.ModelData()

def addImage(index, grid, xrange):
	col = index
	row = 0

	filename = mp1_data.getRadioDirname(5, 2, 25, angles[index], 5) + "/emeasure_ff.fits"

	print filename

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

	if col != 0:
		ax.tick_params(labelleft='off')

	# Add column and row labels
	if row == 0:
		ax.text(0.5, 1.24, fmt.latexify("\\theta_\\mathrm{i} = " + str(angles[col]) + "^\\circ"),
				fontsize=10, horizontalalignment='center', verticalalignment='bottom',
				rotation='horizontal', transform=ax.transAxes)

hrats = [1, 1, 1, 1, 1]
vrats = [1]
fancy_grid = torch.FancyAxesGrid(hrats, vrats, border=(0.02, 0.015, 0.08, 0.16),
								 fig_width=plot_size, dpi=DPI, fig_format=figformat,
								 csize=0.04, cspace=0.012, cpad=0.0175,
								 hspace=-0.000, vspace=0.022)

xrange = [[120], [120], [120], [120], [120]]

for i in range(5):
	addImage(i, fancy_grid, xrange)

fancy_grid.update_cbar()

fancy_grid.save_plot("emeasure-angle.png")