import numpy as np
import matplotlib.ticker as ticker
import warnings
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
plot_size = 5
fontsize = 16
outputfile_qlyc = "mass-vs-mdot.png"

torch.set_font_sizes(fontsize)

data = np.genfromtxt("refdata/zams.txt", skip_header=1)

### Data set up.

mass = data[:,0]
mdot = data[:,8]*1.5768e-26

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
grid = plotter.axes1D((1,1), aspect_ratio=0.75)
grid[0].set_xlabel(plotter.format_label(torch.VarType('M_\star', units='M_{\odot}')))
grid[0].set_ylabel(plotter.format_label(torch.VarType('\dot{M}', units='M_{\odot}\, yr^{-1}', isLog10=True)))

### Plot.
formatter = ticker.ScalarFormatter(useOffset=True, useMathText=True)
formatter.set_powerlimits((0, 1))

grid[0].plot(mass, np.log10(mdot), color='r', linewidth=1)
#grid[0].yaxis.set_major_formatter(formatter)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile_qlyc)

print sys.argv[0] + ': plotted in ' + outputfile_qlyc