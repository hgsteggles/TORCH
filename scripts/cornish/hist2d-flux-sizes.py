import sys
import warnings
import numpy as np
from scipy import ndimage, stats
import matplotlib.pyplot as plt
import math

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
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

cornish_data = cdat.CornishData(args.iden)
star_data = cornish_data.star_data
nstars = len(star_data[:,0])

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
outputfile = cornish_data.dirname + '/hist2d-flux-size' + benstr + '.' + figformat

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotter.ticklength *= 0.5

###	Axes.
asp_rat = 0.6
grid = plotter.axes1D((2,2), aspect_ratio=asp_rat, axes_pad=(0.0, 0.2))

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

### Data
simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + benstr + ".txt", skip_header=1)
cornish_survey = np.genfromtxt("data/cornish/cornish-uchii-distances.txt", skip_header=1, delimiter=',')

phy_sizes = simulated_survey[:,7] * simulated_survey[:,4] * 1000.0 * math.pi / (60.0 * 60.0 * 180.0)
cornish_phy_sizes = cornish_survey[:,11] * cornish_survey[:,33] * 1000.0 * math.pi / (60.0 * 60.0 * 180.0)
############

for iax in range(4):
	grid[iax].set_xlabel(plotter.format_label(torch.VarType('\mathrm{Flux}', units='mJy')))
	grid[iax].set_xlim([2, 10000])
	grid[iax].set_xscale("log")

for iax in range(2):
	grid[iax].set_ylabel(plotter.format_label(torch.VarType('\mathrm{Angular\ Size}', units='\prime\prime')))
	grid[iax].set_ylim([0.0, 15.0])
	grid[iax].set_yticks(np.arange(0, 15, 2))

for iax in range(2, 4):
	grid[iax].set_ylabel(plotter.format_label(torch.VarType('\mathrm{Physical\ Size}', units='pc')))
	grid[iax].set_ylim([0.0, 1.0])


kernal_density_estimate(grid[0], simulated_survey[:,9], simulated_survey[:,7],
						2.0, 10000.0, 0.0, 15.0)
grid[0].scatter(simulated_survey[:,9], simulated_survey[:,7], marker='o',
				  color='k', s=2)

kernal_density_estimate(grid[2], simulated_survey[:,9], phy_sizes,
						2.0, 10000.0, 0.0, 1.0)
grid[2].scatter(simulated_survey[:,9], phy_sizes, marker='o',
				  color='k', s=2)

kernal_density_estimate(grid[1], cornish_survey[:,9], cornish_survey[:,11],
						2.0, 10000.0, 0.0, 15.0)
grid[1].scatter(cornish_survey[:,9], cornish_survey[:,11], marker='o',
				  color='k', s=2)

kernal_density_estimate(grid[3], cornish_survey[:,9], cornish_phy_sizes,
						2.0, 10000.0, 0.0, 1.0)
grid[3].scatter(cornish_survey[:,9], cornish_phy_sizes, marker='o',
				  color='k', s=2)

for iax in range(0, 4, 2):
	grid[iax].text(0.04, 0.94, "Simulated",
							fontsize=fontsize, horizontalalignment='left',
							verticalalignment='top', rotation='horizontal',
							transform=grid[iax].transAxes)

for iax in range(1, 4, 2):
	grid[iax].text(0.04, 0.94, "CORNISH",
							fontsize=fontsize, horizontalalignment='left',
							verticalalignment='top', rotation='horizontal',
							transform=grid[iax].transAxes)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile
