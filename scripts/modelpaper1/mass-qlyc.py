import numpy as np
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
outputfile_qlyc = "mass-vs-qlyc.png"
torch.set_font_sizes(fontsize)

data = np.genfromtxt("refdata/zams.txt", skip_header=1)

### Data set up.

mass = data[:,0]
qlyc = data[:,6]

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
grid = plotter.axes1D((1,1), aspect_ratio=0.75)
grid[0].set_xlabel(plotter.format_label(torch.VarType('M_\star', units='M_{\odot}')))
grid[0].set_ylabel(plotter.format_label(torch.VarType('Q_\\mathrm{Lyc}', units='s^{-1}', isLog10=True)))

### Plot.
grid[0].plot(mass, qlyc, color='r', linewidth=1)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile_qlyc)

print sys.argv[0] + ': plotted in ' + outputfile_qlyc