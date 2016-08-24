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
mp1_data = mp1.ModelData()

DPI = 300
figformat = 'png'
plot_size = 10
torch.set_font_sizes(fontsize=8)
cmap = hgspy.get_grey_cmap(isReversed=True)

snapshots = [10, 20, 30, 40, 50]

def addImage(col, row, grid, xrange):
	filename = mp1_data.getDataFilename(row, 2, snapshots[col])

	ax = grid.grid[col][row]
	cbar_ax = grid.cgrid[col][row]

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

	if col == 0:
		ax.text(-0.15, 0.5, fmt.latexify("M = " + fmt.fmt_mass(mp1_data.masses[row]) + "\ \\mathrm{M_\\odot}"),
				fontsize=10, horizontalalignment='right', verticalalignment='center',
				rotation='vertical', transform=ax.transAxes)
	if row == 0:
		ax.text(0.5, 2.4, fmt.latexify("t = " + fmt.fmt_nolatex(mp1_data.times[col], 0) + "\ \\mathrm{yrs}"),
				fontsize=10, horizontalalignment='center', verticalalignment='bottom',
				rotation='horizontal', transform=ax.transAxes)


hrats = [1, 1, 1, 1, 1]
vrats = [1, 1, 1, 1, 1, 1, 1, 1, 1]
fancy_grid = torch.FancyAxesGrid(hrats, vrats, border=(0.05, 0.01, 0.03, 0.02),
								 fig_width=plot_size, dpi=DPI, fig_format=figformat,
								 csize=0.005, cspace=0.002, cpad=0.03,
								 hspace=0.026, vspace=0.025)

xrange = []
xrange.append([0.040, 0.08, 0.18, 0.28, 0.40, 0.52, 0.6, 0.9, 2.2])
xrange.append([0.035, 0.10, 0.22, 0.34, 0.52, 0.70, 0.9, 1.9, 2.5])
xrange.append([0.038, 0.10, 0.24, 0.40, 0.65, 0.90, 1.0, 2.2, 3.0])
xrange.append([0.040, 0.10, 0.25, 0.48, 0.74, 1.00, 1.3, 2.3, 3.2])
xrange.append([0.040, 0.10, 0.25, 0.54, 0.80, 1.00, 1.6, 2.6, 3.6])

for i in range(5):
	for j in range(9):
		addImage(i, j, fancy_grid, xrange)

fancy_grid.update_cbar()

fancy_grid.save_plot("mass-time-tem-samesize.png")