import sys
import warnings
import numpy as np

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
outputfile = cornish_data.dirname + '/hist-age' + benstr + '.' + figformat

### Data
simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + benstr + ".txt", skip_header=1)

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotter.ticklength *= 0.5

###	Axes.
asp_rat = 1.0
grid = plotter.axes1D((1,1), aspect_ratio=asp_rat)

grid[0].set_xlabel(plotter.format_label(torch.VarType('\mathrm{Age}', units='kyr')))

grid[0].set_ylabel(plotter.format_label(torch.VarType('N')))

grid[0].set_xlim([0, 200])

bins = np.arange(0, 200, 5)
kx = dict(linewidth=1.5, label="Simulated", color='r', linestyle='--')
plotter.histstep(grid[0], simulated_survey[:,3], bins, **kx)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile
