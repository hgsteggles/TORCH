import sys
import warnings
import numpy as np

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Plots location histograms of galaxy population.')
parser.add_argument('filename', metavar='filename', type=str, help='starpop file.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

DPI = 300
figformat = 'png'
plot_size = 10.0
fontsize = 16
torch.set_font_sizes(fontsize)

outputfile = 'hist-locations-full' + '.' + figformat

### Data
simulated_survey = np.genfromtxt(args.filename, skip_header=1)
cornish_survey = np.genfromtxt("data/cornish/cornish-uchiis.txt", skip_header=1, delimiter=',')
cornish_distances = np.genfromtxt("data/cornish/cornish-distances.txt", skip_header=1)

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotter.ticklength *= 0.5

###	Axes.
asp_rat = 1.0
grid = plotter.axes1D((1,3), aspect_ratio=asp_rat)

grid[0].set_xlabel(plotter.format_label(torch.VarType('l', units='deg')))
grid[1].set_xlabel(plotter.format_label(torch.VarType('b', units='deg')))
grid[2].set_xlabel(plotter.format_label(torch.VarType('d', units='kpc')))

grid[0].set_ylabel(plotter.format_label(torch.VarType('N')))

grid[0].set_xlim([10.0, 65.0])
grid[1].set_xlim([-1, 1])
grid[2].set_xlim([0.0, 20.0])

kx1 = dict(linewidth=1.5, label="CORNISH")
kx2 = dict(linewidth=1.5, label="Simulated", linestyle='--')
plotter.histstep(grid[0], cornish_survey[:,1], np.arange(10, 65, 2), **kx1)
plotter.histstep(grid[0], simulated_survey[:,15], np.arange(10, 65, 2), **kx2)
plotter.histstep(grid[1], cornish_survey[:,2], np.arange(-1, 1, 0.1), **kx1)

plotter.histstep(grid[1], simulated_survey[:,16], np.arange(-1, 1, 0.1), **kx2)
plotter.histstep(grid[2], cornish_distances[:,1], np.arange(0, 20, 1), **kx1)
plotter.histstep(grid[2], simulated_survey[:,4], np.arange(0, 20, 1), **kx2)

### Legend
handles, labels = grid[1].get_legend_handles_labels()
legend = grid[1].legend(handles, labels, loc=1)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile