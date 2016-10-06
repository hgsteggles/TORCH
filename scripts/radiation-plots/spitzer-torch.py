import warnings

import numpy as np
import sys
import math

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

outputfile = 'spitzer-torch' + '.' + figformat

inputfile = "data/spitzer-torch/if.csv"


###	Data set up.
S2YR = 1.0 / (365.0 * 24.0 * 60.0 * 60.0)
data = np.genfromtxt(inputfile, delimiter=',')

t = data[:,0] #[yrs]
R = data[:,1] * 3.09e18 #[cm]

nh = 400.0
alphaB = 2.59e-13
Q = 1.0e49
ts = 3.0e12 * S2YR
RS = np.power(3.0 * Q / (4.0 * math.pi * nh * nh * alphaB), 1.0 / 3.0)
R_spit = RS * np.power(np.ones(t.shape) + (7.0 / 4.0) * t / ts, 4.0 / 7.0)
R_raga = RS * np.power(np.ones(t.shape) + (7.0 / 4.0) * np.sqrt(4.0 / 3.0) * t / ts, 4.0 / 7.0)

error_spit = np.abs((R_spit - R) / R_spit)
error_raga = np.abs((R_raga - R) / R_raga)

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
asp_rat = 7.0 / 16.0
grid = plotter.axes1D((2,1), aspect_ratio=asp_rat)
grid[0].yaxis.set_ticks(np.arange(0.0, 0.11, 0.02))
grid[1].yaxis.set_ticks(np.arange(0.0, 8.0, 1.0))
grid[0].set_xlim([0.0, 16.0])
grid[0].set_ylim([0.0, 0.10])
grid[1].set_ylim([0, 8])
grid[1].set_xlabel(plotter.format_label(torch.VarType('t\ /\ t_\\mathrm{s}')))
grid[0].set_ylabel(plotter.format_label(torch.VarType('\\mathrm{Relative\ Error}')))
grid[1].set_ylabel(plotter.format_label(torch.VarType('R_\mathrm{IF}\ /\ R_\mathrm{st}')))

### Plot.
kx = dict(linewidth=1.5)

grid[1].plot(t / ts, R / RS, color='black', label='Simulation', **kx)
grid[1].plot(t / ts, R_spit / RS, color='blue', label='Spitzer', **kx)
grid[1].plot(t / ts, R_raga / RS, color='red', label='Raga', **kx)

grid[0].plot(t / ts, error_spit, color='blue', label='Spitzer', **kx)
grid[0].plot(t / ts, error_raga, color='red', label='Raga', **kx)

### Legend
handles, labels = grid[1].get_legend_handles_labels()
legend = grid[0].legend(handles, labels)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted ' + outputfile

