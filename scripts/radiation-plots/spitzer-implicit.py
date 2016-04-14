import warnings

import numpy as np
import sys

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

outputfile = 'spitzer' + '.' + figformat

inputfile = "data/spitzer-implicit/IF_128.dat"


###	Data set up.
data = np.genfromtxt(inputfile)

t = data[:,0]
RS = data[:,1]
RS_a = data[:,2]
error = data[:,3]

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
asp_rat = 7.0 / 16.0
grid = plotter.axes1D((2,1), aspect_ratio=asp_rat)
grid[0].yaxis.set_ticks(np.arange(0.0, 0.06, 0.01))
grid[1].yaxis.set_ticks(np.arange(0.0, 7.0, 1.0))
grid[0].set_xlim([0.0, 16.0])
grid[0].set_ylim([0.0, 0.05])
grid[1].set_ylim([0, 7])
grid[1].set_xlabel(plotter.format_label(torch.VarType('t\ /\ t_\\mathrm{s}', False)))
grid[0].set_ylabel(plotter.format_label(torch.VarType('\\mathrm{Relative\ Error}', False)))
grid[1].set_ylabel(plotter.format_label(torch.VarType('R_\mathrm{IF}\ /\ R_\mathrm{s}', False)))

### Plot.
kx = dict(linewidth=1.5)

grid[1].plot(t, RS, color='red', label='Implicit', **kx)
grid[1].plot(t, RS_a, color='black', label='Analytical', **kx)

grid[0].plot(t, error, color='red', label='Implicit', **kx)

### Legend
handles, labels = grid[1].get_legend_handles_labels()
legend = grid[0].legend(handles, labels)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted multi-shadow in ' + outputfile

