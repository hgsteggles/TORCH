import sys
import warnings
import numpy as np

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
plot_size = 5.0
fontsize = 10
torch.set_font_sizes(fontsize)

benstr = ""
if args.ben:
	benstr = "-ben"
if args.harry:
	benstr = "-harry"
outputfile = 'multi-hist-flux' + benstr + '.' + figformat

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotter.ticklength *= 0.5

###	Axes.
asp_rat = 0.4
grid = plotter.axes1D((5,1), aspect_ratio=asp_rat)

for irow in range(5):
	cornish_data = cdat.CornishData(irow)
	star_data = cornish_data.star_data
	nstars = len(star_data[:,0])

	### Data
	simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + benstr + ".txt", skip_header=1)
	cornish_survey = np.genfromtxt("data/cornish/cornish-uchiis.txt", skip_header=1, delimiter=',')

	grid[irow].set_xlabel(plotter.format_label(torch.VarType('\mathrm{Flux}', units='mJy')))
	grid[irow].set_ylabel(plotter.format_label(torch.VarType('P')))
	grid[irow].set_xlim([2, 20000])
	grid[irow].set_ylim([0, 0.2])
	grid[irow].set_xscale("log")

	if irow == 4:
		grid[irow].set_yticks(np.arange(0, 0.21, 0.05))
	else:
		grid[irow].set_yticks(np.arange(0.05, 0.21, 0.05))
	#grid[irow].set_yticks(range(0, 45, 5))

	grid[irow].text(0.04, 0.94, fmt.latexify("n_\\star = " + fmt.fmt_power(mp1_data.densities[irow], '{:3.1f}', 4) + "\ \\mathrm{cm^{-3}}"),
				fontsize=fontsize, horizontalalignment='left', verticalalignment='top',
				rotation='horizontal', transform=grid[irow].transAxes)

	kx1 = dict(linewidth=1.5, label="CORNISH", color='b')
	kx2 = dict(linewidth=1.5, label="Simulated", color='r', linestyle='--')

	bins = np.logspace(0.3, 4.3, 22)
	bincentres = np.sqrt(bins[:-1] * bins[1:])

	plotter.histstep(grid[irow], cornish_survey[:,9], bins, errorcentres=bincentres, normed=True, **kx1)
	plotter.histstep(grid[irow], simulated_survey[:,9], bins, normed=True, **kx2)

### Legend
handles, labels = grid[0].get_legend_handles_labels()
legend = grid[0].legend(handles, labels, loc=1)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile
