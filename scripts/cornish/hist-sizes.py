import sys
import warnings
import numpy as np
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
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

load_src("cdat", "../torchpack/cornishdata.py")
import cdat
cornish_data = cdat.CornishData(args.iden)
star_data = cornish_data.star_data
nstars = len(star_data[:,0])

DPI = 300
figformat = 'png'
plot_size = 10.0
fontsize = 16
torch.set_font_sizes(fontsize)

benstr = ""
if args.ben:
	benstr = "-ben"
if args.harry:
	benstr = "-harry"
outputfile = cornish_data.dirname + '/hist-sizes' + benstr + '.' + figformat

### Data
simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + benstr + ".txt", skip_header=1)
cornish_survey = np.genfromtxt("data/cornish/cornish-uchii-distances.txt", skip_header=1, delimiter=',')

phy_sizes = simulated_survey[:,7] * simulated_survey[:,4] * 1000.0 * math.pi / (60.0 * 60.0 * 180.0)
cornish_phy_sizes = cornish_survey[:,11] * cornish_survey[:,33] * 1000.0 * math.pi / (60.0 * 60.0 * 180.0)

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotter.ticklength *= 0.5

###	Axes.
asp_rat = 1.0
grid = plotter.axes1D((1,2), aspect_ratio=asp_rat)

grid[0].set_xlabel(plotter.format_label(torch.VarType('\mathrm{Angular\ Size}', units='arcsec')))
grid[1].set_xlabel(plotter.format_label(torch.VarType('\mathrm{Physical\ Size}', units='pc')))

grid[0].set_ylabel(plotter.format_label(torch.VarType('N')))

grid[0].set_xlim([0.0, 24.0])
grid[1].set_xlim([0.0, 1.0])

### Plot.
#grid[0].hist(simulated_survey[:,7], bins=np.arange(0, 24, 1), label="Model")
#grid[1].hist(phy_sizes, bins=np.arange(0, 1, 0.04), label="Model")

#grid[0].hist(cornish_survey[:,11], bins=np.arange(0, 24, 1), alpha=0.5, color='r', label="CORNISH")
#grid[1].hist(cornish_phy_sizes, bins=np.arange(0, 1, 0.04), alpha=0.5, color='r', label="CORNISH")

kx1 = dict(linewidth=1.5, label="CORNISH", linestyle='-')
kx2 = dict(linewidth=1.5, label="Simulated", linestyle='--')
plotter.histstep(grid[0], cornish_survey[:,11], np.arange(0, 24, 1), **kx1)
plotter.histstep(grid[0], simulated_survey[:,7], np.arange(0, 24, 1), **kx2)
plotter.histstep(grid[1], cornish_phy_sizes, np.arange(0, 1, 0.04), **kx1)
plotter.histstep(grid[1], phy_sizes, np.arange(0, 1, 0.04), **kx2)

### Legend
handles, labels = grid[1].get_legend_handles_labels()
legend = grid[1].legend(handles, labels, loc=1)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile
