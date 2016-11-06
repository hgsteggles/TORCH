import sys
import warnings
import numpy as np

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

DPI = 300
figformat = 'png'
plot_size = 5.0
fontsize = 12
torch.set_font_sizes(fontsize)
dir = "data/stromgren-compare/"

outputfile = 'plot-stromgren' + '.' + figformat

inputfile = []
inputfile.append(dir + "front_I_256x1x1_ts100.txt")
inputfile.append(dir + "front_E1_256x1x1_ts100.txt")
inputfile.append(dir + "front_E2_256x1x1_ts100.txt")


###	Data set up.
data = []
for i in range(len(inputfile)):
	data.append(np.genfromtxt(inputfile[i]))

IF_a = data[0][:,2]

IF = []
R = []
E = []
for i in range(len(data)):
	IF.append(data[i][:,1])
	R.append(data[i][:,0])
	E.append(data[i][:,3] - 1)

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
asp_rat = 7.0 / 16.0
grid = plotter.axes1D((2,1), aspect_ratio=asp_rat)
grid[1].yaxis.set_ticks(np.arange(0.0, 1.8, 0.2))
grid[0].yaxis.set_ticks(np.arange(-0.2, 1.2, 0.2))
grid[0].set_xlim([0.0, 3.8])
grid[0].set_ylim([-0.19, 1.0])
grid[1].set_ylim([0, 1.6])
grid[1].set_xlabel(plotter.format_label(torch.VarType('t\ /\ t_\\mathrm{rec}', isLog10=False)))
grid[0].set_ylabel(plotter.format_label(torch.VarType('\mathrm{Relative\ Error}', isLog10=False)))
grid[1].set_ylabel(plotter.format_label(torch.VarType('R_\mathrm{IF}\ /\ \mathrm{pc}', isLog10=False)))

### Plot.
kx = dict(linewidth=1.0)
kx2 = dict(linewidth=0.8)

grid[1].plot(R[0], IF[0], color='green', label='Implicit', linestyle='-', **kx2)
grid[1].plot(R[1], IF[1], color='red', label='Explicit 1st Order', linestyle='-.', **kx)
grid[1].plot(R[2], IF[2], color='blue', label='Explicit 2nd Order', linestyle='-', **kx)
grid[1].plot(R[0], IF_a, color='black', label='Analytical', linestyle='--', **kx)

grid[0].plot(R[0], E[0], color='green', label='Implicit', linestyle='-', **kx2)
grid[0].plot(R[1], E[1], color='red', label='Explicit 1st Order', linestyle='-.', **kx)
grid[0].plot(R[2], E[2], color='blue', label='Explicit 2nd Order', linestyle='-', **kx)

### Legend
handles, labels = grid[1].get_legend_handles_labels()
legend = grid[0].legend(handles, labels, fontsize=10)

legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile

