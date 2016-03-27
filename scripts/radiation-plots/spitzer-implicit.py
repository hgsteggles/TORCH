import warnings

import numpy as np
import sys
sys.path.insert(0, '/home/harry/Documents/Projects/PhD/Torch/scripts/blah')
import torchpack.torch as torch
import torchpack.hgspy as hgspy

DPI = 300
figformat = 'png'
plot_size = 2.5

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
grid = plotter.axes1D((2,1))
grid[0].yaxis.set_ticks(np.arange(0.0, 0.06, 0.01))
grid[1].yaxis.set_ticks(np.arange(0.0, 7.0, 1.0))
grid[0].set_xlim([0.0, 16.0])
grid[0].set_ylim([0.0, 0.05])
grid[1].set_ylim([0, 7])
grid[1].set_xlabel(plotter.format_label(torch.VarType('t\ /\ t_s', False)))
grid[0].set_ylabel(plotter.format_label(torch.VarType('Relative\ Error', False)))
grid[1].set_ylabel(plotter.format_label(torch.VarType('R_{IF}\ /\ R_s', False)))

### Plot.
kx = dict(linewidth=0.5)

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

