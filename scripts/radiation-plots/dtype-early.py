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
plot_size = 5
fontsize = 12
torch.set_font_sizes(fontsize)

outputfile = 'dtype-1d-early-2000' + '.' + figformat

inputfile = "data/D-type/early-1D-hrhr/if.csv"


###	Data set up.
S2YR = 1.0 / (365.0 * 24.0 * 60.0 * 60.0)
PC2CM = 3.09e18
RGAS = 8.3144598e7

data = np.genfromtxt(inputfile, delimiter=',')

t = data[:,0] / S2YR #[yrs]
R = data[:,1] * PC2CM #[cm]

nh = 3112.3
alphaB = 2.59e-13
Q = 1.0e49
RS = np.power(3.0 * Q / (4.0 * math.pi * nh * nh * alphaB), 1.0 / 3.0)

def soundSpeed(T, mu):
	return math.sqrt(RGAS * T / mu)

ci = soundSpeed(10000.0, 0.5)
co = soundSpeed(100.0, 1.0)

cico = ci / co
ts = RS / ci
R_stag = (4.0 / 3.0)**(2.0 / 3.0) * (cico)**(4.0 / 3.0) * RS
t_stag = (4.0 / 7.0) * ((3.0 / 4.0)**0.5) * (((4.0 / 3.0)**(7.0 / 6.0)) * (cico**(7.0 / 3.0)) - 1.0) * ts

R_spit = RS * np.power(np.ones(t.shape) + (7.0 / 4.0) * t / ts, 4.0 / 7.0)
R_raga = RS * np.power(np.ones(t.shape) + (7.0 / 4.0) * np.sqrt(4.0 / 3.0) * t / ts, 4.0 / 7.0)

error_spit = np.abs((R_spit - R) / R_spit)
error_raga = np.abs((R_raga - R) / R_raga)

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
asp_rat = 7.0 / 16.0
grid = plotter.axes1D((2,1), aspect_ratio=asp_rat)
grid[0].yaxis.set_ticks(np.arange(0.01, 0.08, 0.01))
grid[1].yaxis.set_ticks(np.arange(0.0, 1.6, 0.2))
grid[0].set_xlim([0.0, 0.14])
grid[1].set_ylim([0, 1.4])
grid[1].set_xlabel(plotter.format_label(torch.VarType('t', units='Myr')))
grid[0].set_ylabel(plotter.format_label(torch.VarType('\\mathrm{Relative\ Error}')))
grid[1].set_ylabel(plotter.format_label(torch.VarType('R_\mathrm{IF}', units='pc')))

### Plot.
kx = dict(linewidth=1.2)
kx2 = dict(linewidth=0.8)

grid[1].plot(t * S2YR / 1.0e6, R / PC2CM, color='black', label='Simulation', linestyle='-', **kx2)
grid[1].plot(t * S2YR / 1.0e6, R_spit / PC2CM, color='blue', label='Spitzer', linestyle='--', **kx)
grid[1].plot(t * S2YR / 1.0e6, R_raga / PC2CM, color='red', label='Hosokawa-Inutsuka', linestyle='-.', **kx)

grid[0].plot(t * S2YR / 1.0e6, error_spit, color='blue', label='Spitzer', linestyle='--', **kx)
grid[0].plot(t * S2YR / 1.0e6, error_raga, color='red', label='Hosokawa-Inutsuka', linestyle='-.', **kx)

### Legend
handles, labels = grid[1].get_legend_handles_labels()
legend = grid[0].legend(handles, labels, loc='right', fontsize=10)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted ' + outputfile

