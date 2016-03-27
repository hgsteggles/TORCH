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
cmap = hgspy.get_temperature_cmap()

fileprefix = "data/model_paper1/set2/"

offgrid = np.fromfile(fileprefix + "offgrid_times", dtype=int, sep='\n')

masses = [6, 9, 12, 15, 20, 30, 40, 70, 120]
densities = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]

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


def getFilename(i):
	padi = "%02d" % ((i+1),)
	return fileprefix + "data_" + padi + "/data2D_025.txt"

def addImage(index, grid, xrange):
	col = int(index / 9)
	row = index % 9

	filename = getFilename(index)

	ax = grid.grid[col][row]
	cbar_ax = grid.cgrid[col][row]


	if offgrid[index] > 25:
		data = torch.CFD_Data(filename, axial=True)

		vs = data.get_var(torch.VarType("tem", True))
		vsi = data.interpolate(vs, "linear")
		vsmin = vs.min()
		vsmax = vs.max()

		im = torch.image(data.x[0], data.x[1],
							 vsi, vsmin, vsmax, "linear", cmap,
							 ax)

		#cb = cbar_ax.colorbar(im, format=ticker.FuncFormatter(fmt))
		formatter = ticker.ScalarFormatter(useOffset=True, useMathText=True)
		formatter.set_powerlimits((0, 1))

		cb = grid.fig.colorbar(im, cax=cbar_ax, format=ticker.FuncFormatter(fmt), orientation='horizontal')
		#cb.ax.xaxis.set_ticks(np.arange(0, 1.01, 0.5))
		#labels = cb.ax.get_xticklabels()
		#labels[0] = ""
		#labels[2] = ""
		#cb.ax.set_xticklabels(labels)
	else:
		grid.set_visible(col, row, False)

	if grid.is_visible(col, row):
		zcentre = (110.0 / 200.0) * data.dx * data.ny

		zmin = zcentre - xrange[col][row]/2.0
		zmax = zcentre + xrange[col][row]/2.0

		ax.set_xlim([-xrange[col][row]/2.0, xrange[col][row]/2.0])
		ax.set_ylim([zmin, zmax])
		ax.xaxis.set_ticks(np.arange(-xrange[col][row]/2.0, (xrange[col][row]/2.0) + 0.0001, xrange[col][row]/4.0))
		ax.yaxis.set_ticks(np.arange(zmin, zmax + 0.0001, xrange[col][row]/4.0))

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
		ax.text(-0.22, 0.5, latexify("M = " + fmt_mass(masses[row]) + "\ \\mathrm{M_\\odot}"),
				fontsize=10, horizontalalignment='right', verticalalignment='center',
				rotation='vertical', transform=ax.transAxes)
	if row == 0:
		ax.text(0.5, 1.24, latexify("n_\\star = " + fmt_nolatex(densities[col], 0) + "\ \\mathrm{cm^{-3}}"),
				fontsize=10, horizontalalignment='center', verticalalignment='bottom',
				rotation='horizontal', transform=ax.transAxes)


hrats = [1, 1, 1, 1, 1]
vrats = [1, 1, 1, 1, 1, 1, 1, 1]
fancy_grid = torch.FancyAxesGrid(hrats, vrats, border=(0.05, 0.01, 0.03, 0.02),
								 fig_width=plot_size, dpi=DPI, fig_format=figformat,
								 csize=0.005, cspace=0.002, cpad=0.03,
								 hspace=0.026, vspace=0.025)

xrange = []
xrange.append([0.10, 0.24, 0.46, 0.70, 0.90, 1.0, 1, 1])
xrange.append([0.06, 0.15, 0.38, 0.60, 0.80, 0.9, 1, 1])
xrange.append([0.04, 0.11, 0.25, 0.45, 0.65, 0.9, 1, 1])
xrange.append([0.03, 0.06, 0.20, 0.35, 0.55, 0.8, 1, 1])
xrange.append([0.02, 0.05, 0.14, 0.25, 0.45, 0.6, 0.8, 1])

for i in range(45):
	if (i + 1) % 9 != 0:
		addImage(i, fancy_grid, xrange)

fancy_grid.update_cbar()

fancy_grid.save_plot("mass-density-tem-samesize.png")