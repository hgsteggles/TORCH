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
cmap = hgspy.get_temperature_cmap()

def addImage(index, grid, xrange):
	col = int(index / 9)
	row = index % 9

	filename = mp1_data.getDataFilename(row, col, 25)

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
	def fmt_mass(x):
		a = '{:.0f}'.format(x)
		return '{}'.format(a)

	if col == 0:
		ax.text(-0.22, 0.5, fmt.latexify("M = " + fmt.fmt_mass(mp1_data.masses[row]) + "\ \\mathrm{M_\\odot}"),
				fontsize=10, horizontalalignment='right', verticalalignment='center',
				rotation='vertical', transform=ax.transAxes)
	if row == 0:
		ax.text(0.5, 1.24, fmt.latexify("n_\\star = " + fmt.fmt_nolatex(mp1_data.densities[col], 0) + "\ \\mathrm{cm^{-3}}"),
				fontsize=10, horizontalalignment='center', verticalalignment='bottom',
				rotation='horizontal', transform=ax.transAxes)


hrats = [1, 1, 1, 1, 1]
vrats = [1, 1, 1, 1, 1, 1, 1, 1, 1]
fancy_grid = torch.FancyAxesGrid(hrats, vrats, border=(0.05, 0.01, 0.03, 0.02),
								 fig_width=plot_size, dpi=DPI, fig_format=figformat,
								 csize=0.005, cspace=0.002, cpad=0.03,
								 hspace=0.026, vspace=0.025)

xrange = []
xrange.append([15, 32, 62, 90, 134, 190, 260, 380, 800])
xrange.append([10, 22, 48, 80, 110, 140, 190, 300, 500])
xrange.append([6, 15, 34, 60, 90, 120, 150, 280, 400])
xrange.append([4, 9, 24, 46, 66, 100, 120, 180, 250])
xrange.append([2.52, 6, 16, 30, 54, 80, 100, 150, 190])

for i in range(5):
	for j in range(9):
		xrange[i][j] = xrange[i][j] * 0.1 / 15.0

#xrange = []
#xrange.append([0.10, 0.24, 0.46, 0.70, 0.90, 1.0, 1, 1])
#xrange.append([0.06, 0.15, 0.38, 0.60, 0.80, 0.9, 1, 1])
#xrange.append([0.04, 0.11, 0.25, 0.45, 0.65, 0.9, 1, 1])
#xrange.append([0.03, 0.06, 0.20, 0.35, 0.55, 0.8, 1, 1])
#xrange.append([0.02, 0.05, 0.14, 0.25, 0.45, 0.6, 0.8, 1])

for i in range(45):
	addImage(i, fancy_grid, xrange)

fancy_grid.update_cbar()

fancy_grid.save_plot("mass-density-tem-samesize.png")