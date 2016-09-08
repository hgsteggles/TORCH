import sys
import warnings
import numpy as np
from scipy import ndimage, stats
import matplotlib.pyplot as plt

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('-ben', action='store_true', help='Use simple ben stromgren survey.')
parser.add_argument('-harry', action='store_true', help='Use simple harry stromgren survey.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
import torch
load_src("hgspy", "../torchpack/hgspy.py")
import hgspy
load_src("cdat", "../torchpack/cornishdata.py")
import cdat
load_src("fmt", "../torchpack/formatting.py")
import fmt
load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData()

DPI = 300
figformat = 'png'
plot_size = 4.0
fontsize = 13
torch.set_font_sizes(fontsize)
colormap = plt.cm.bone_r

benstr = ""
if args.ben:
	benstr = "-ben"
if args.harry:
	benstr = "-harry"
outputfile = 'multi-hist2d-flux-size' + benstr + '.' + figformat

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotter.ticklength *= 0.5

###	Axes.
asp_rat = 0.6
grid = plotter.axes1D((5,2), aspect_ratio=asp_rat)

def kernal_density_estimate(ax, pos_x, pos_y, xmin, xmax, ymin, ymax):
	xmin = np.log10(xmin)
	xmax = np.log10(xmax)
	pos_x = np.log10(pos_x)
	#print pos_x
	X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
	positions = np.vstack([X.ravel(), Y.ravel()])
	values = np.vstack([pos_x, pos_y])
	kernel = stats.gaussian_kde(values)
	Z = np.reshape(kernel(positions).T, X.shape)
	X = np.power(10.0, X)
	ax.contourf(X, Y, Z, cmap=colormap, alpha=0.8)


for irow in range(5):
	cornish_data = cdat.CornishData(irow)
	star_data = cornish_data.star_data
	nstars = len(star_data[:,0])

	### Data
	simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + benstr + ".txt", skip_header=1)
	cornish_survey = np.genfromtxt("data/cornish/cornish-uchiis.txt", skip_header=1, delimiter=',')
	############

	for icol in range(2):
		grid[2 * irow + icol].set_xlabel(plotter.format_label(torch.VarType('\mathrm{Flux}', units='mJy')))
		grid[2 * irow + icol].set_ylabel(plotter.format_label(torch.VarType('\mathrm{Angular\ Size}', units='\prime\prime')))
		grid[2 * irow + icol].set_xlim([2, 10000])
		grid[2 * irow + icol].set_ylim([0.0, 15.0])

		if irow == 4:
			for icol in range(2):
				grid[2 * irow + icol].set_yticks(np.arange(0, 15, 2))
		else:
			for icol in range(2):
				grid[2 * irow + icol].set_yticks(np.arange(2, 15, 2))

	grid[2 * irow + 0].text(0.04, 0.94, fmt.latexify("n_\\star = " + fmt.fmt_power(mp1_data.densities[irow], '{:3.1f}', 4) + "\ \\mathrm{cm^{-3}}"),
					fontsize=fontsize, horizontalalignment='left', verticalalignment='top',
					rotation='horizontal', transform=grid[2 * irow + 0].transAxes)

	kernal_density_estimate(grid[2 * irow + 0], simulated_survey[:,9], simulated_survey[:,7], 2.0, 10000.0, 0.0, 15.0)

	grid[2 * irow + 0].scatter(simulated_survey[:,9], simulated_survey[:,7],
							   marker='o', color='k', s=2)

	if irow == 4:
		grid[2 * irow + 1].text(0.04, 0.94, "CORNISH",
					fontsize=fontsize, horizontalalignment='left', verticalalignment='top',
					rotation='horizontal', transform=grid[2 * irow + 1].transAxes)

		kernal_density_estimate(grid[2 * irow + 1], cornish_survey[:,9], cornish_survey[:,11], 2.0, 10000.0, 0.0, 15.0)

		grid[2 * irow + 1].scatter(cornish_survey[:,9], cornish_survey[:,11],
								   marker='o', color='k', s=2)
	else:
		grid[2 * irow + 1].set_visible(False)

	grid[2 * irow + 0].set_xscale("log")
	grid[2 * irow + 1].set_xscale("log")

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile
