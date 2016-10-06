import sys
import warnings
import numpy as np
import math

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

benstr = ""
if args.ben:
	benstr = "-ben"
if args.harry:
	benstr = "-harry"
outputfile = 'multi-hist-sizes' + benstr + '.' + figformat

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotter.ticklength *= 0.5

###	Axes.
asp_rat = 0.6
grid = plotter.axes1D((5,2), aspect_ratio=asp_rat)

normed = False

ymaxes = 5 * [0.4]
nyticks = 5 * [4]

if not normed:
	ymaxes = [100, 100, 100, 160, 300]
	nyticks = [4, 4, 4, 4, 4]

for irow in range(5):
	cornish_data = cdat.CornishData(irow)
	star_data = cornish_data.star_data
	nstars = len(star_data[:,0])

	### Data
	simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + benstr + ".txt", skip_header=1)
	cornish_survey = np.genfromtxt("data/cornish/cornish-uchii-distances.txt", skip_header=1, delimiter=',')

	phy_sizes = simulated_survey[:,7] * simulated_survey[:,4] * 1000.0 * math.pi / (60.0 * 60.0 * 180.0)
	cornish_phy_sizes = cornish_survey[:,11] * cornish_survey[:,33] * 1000.0 * math.pi / (60.0 * 60.0 * 180.0)

	grid[2 * irow + 0].set_xlabel(plotter.format_label(torch.VarType('l', units='deg')))
	grid[2 * irow + 1].set_xlabel(plotter.format_label(torch.VarType('b', units='deg')))

	grid[2 * irow + 0].set_xlabel(plotter.format_label(torch.VarType('\mathrm{Angular\ Size}', units='\prime\prime')))
	grid[2 * irow + 1].set_xlabel(plotter.format_label(torch.VarType('\mathrm{Physical\ Size}', units='pc')))

	ylabel = 'P' if normed else 'N'
	grid[2 * irow + 0].set_ylabel(plotter.format_label(torch.VarType(ylabel)))

	grid[2 * irow + 0].set_xlim([0.0, 24.0])
	grid[2 * irow + 1].set_xlim([0.0, 1.0])

	ymax = ymaxes[irow]
	nyt = nyticks[irow]

	grid[2 * irow + 0].set_ylim([0, ymax])
	grid[2 * irow + 1].set_ylim([0, ymax])

	if irow == 4:
		for icol in range(2):
			grid[2 * irow + icol].set_yticks(np.arange(0, ymax +  (ymax - 1.0e-8)/ float(nyt), ymax / float(nyt)))
	else:
		for icol in range(2):
			grid[2 * irow + icol].set_yticks(np.arange(ymax / float(nyt), ymax +  (ymax - 1.0e-8)/ float(nyt), ymax / float(nyt)))

	kx1 = dict(linewidth=1.5, label="CORNISH", color='b')
	kx2 = dict(linewidth=1.5, label="Simulated", color='r', linestyle='--')

	grid[2 * irow + 0].text(0.04, 0.94, fmt.latexify("n_\\star = " + fmt.fmt_power(mp1_data.densities[irow], '{:3.1f}', 4) + "\ \\mathrm{cm^{-3}}"),
				fontsize=fontsize, horizontalalignment='left', verticalalignment='top',
				rotation='horizontal', transform=grid[2 * irow + 0].transAxes)

	bins0 = np.arange(0, 24, 1)
	bins1 = np.arange(0, 1, 0.04)

	bincentres0 = 0.5 * (bins0[:-1] + bins0[1:])
	bincentres1 = 0.5 * (bins1[:-1] + bins1[1:])

	plotter.histstep(grid[2* irow + 0], cornish_survey[:,11], bins0, errorcentres=bincentres0, normed=normed, **kx1)
	plotter.histstep(grid[2* irow + 0], simulated_survey[:,7], bins0, normed=normed, **kx2)
	plotter.histstep(grid[2* irow + 1], cornish_phy_sizes, bins1, errorcentres=bincentres1, normed=normed, **kx1)
	plotter.histstep(grid[2* irow + 1], phy_sizes, bins1, normed=normed, **kx2)

### Legend
handles, labels = grid[1].get_legend_handles_labels()
legend = grid[1].legend(handles, labels, loc=1)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile
