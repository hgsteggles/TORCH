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

snapshots = [10, 20, 30, 40, 50]
param_index = range(19, 28, 1)

offgrid = np.fromfile(fileprefix + "offgrid_times", dtype=int)

masses = [6, 9, 12, 15, 20, 30, 40, 70, 120]
times = [2.0e4, 4.0e4, 6.0e4, 8e4, 10e4]

def latexify(str):
	return r'${}$'.format(str)

def fmt_nolatex(x, pos):
	if x >= -100 and x <= 100:
		a = '{:.1f}'.format(x)
		return '{}'.format(a)
	else:
		a, b = '{:.1e}'.format(x).split('e')
		b = int(b)
		return '{} \\times 10^{{{}}}'.format(a, b)

def fmt(x, pos):
	return latexify(fmt_nolatex(x, pos))


def getFilename(i, j):
	padi = "%03d" % (snapshots[i],)
	padj = "%02d" % (param_index[j],)

	return fileprefix + "data_" + padj + "/radio_" + padi + "/emeasure_ff.fits"

def addImage(col, row, grid, hrange, vrange):
	filename = getFilename(col, row)

	ax = grid.grid[col][row]
	cbar_ax = grid.cgrid[col][row]

	if offgrid[param_index[row] - 1] > snapshots[col]:
		hdu_list = fits.open(filename)
		image_data = hdu_list[0].data
		image_data = ndimage.gaussian_filter(image_data, sigma=(5, 5), order=0)

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

		cb = grid.fig.colorbar(im, cax=cbar_ax, format=ticker.FuncFormatter(fmt),
							   orientation='horizontal')
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

	if grid.is_visible(col, row):
		ax.set_xlim([-hrange/2, hrange/2])
		ax.set_ylim([-vrange/2, vrange/2])
		ax.xaxis.set_ticks(np.arange(-hrange/2, (hrange/2) + 1, hrange/4))
		ax.yaxis.set_ticks(np.arange(-vrange/2, (vrange/2) + 1, vrange/4))

		ax.xaxis.get_major_ticks()[0].label1.set_visible(False)
		ax.xaxis.get_major_ticks()[2].label1.set_visible(False)
		ax.xaxis.get_major_ticks()[-1].label1.set_visible(False)

		ax.yaxis.get_major_ticks()[0].label1.set_visible(False)
		ax.yaxis.get_major_ticks()[2].label1.set_visible(False)
		ax.yaxis.get_major_ticks()[-1].label1.set_visible(False)

	# Add column and row labels
	def fmt_mass(x):
		a = '{:.0f}'.format(x)
		return '{}'.format(a)

	if col == 0:
		ax.text(-0.15, 0.5, latexify("M = " + fmt_mass(masses[row]) + "\ \\mathrm{M_\\odot}"),
				fontsize=10, horizontalalignment='right', verticalalignment='center',
				rotation='vertical', transform=ax.transAxes)
	if row == 0:
		ax.text(0.5, 2.4, latexify("t = " + fmt_nolatex(times[col], 0) + "\ \\mathrm{yrs}"),
				fontsize=10, horizontalalignment='center', verticalalignment='bottom',
				rotation='horizontal', transform=ax.transAxes)


hrats = [160.0, 160.0, 160.0, 160.0, 160.0]
vrats = [20.0, 40.0, 60.0, 80.0, 120.0, 160.0, 160.0, 160.0]
fancy_grid = torch.FancyAxesGrid(hrats, vrats, border=(0.05, 0.01, 0.03, 0.02),
								 fig_width=plot_size, dpi=DPI, fig_format=figformat,
								 csize=0.005, cspace=0.002, hspace=0.002, vspace=0.025)

for i in range(5):
	for j in range(8):
		addImage(i, j, fancy_grid, hrats[i], vrats[j])
fancy_grid.update()

fancy_grid.save_plot("mass-time.png")