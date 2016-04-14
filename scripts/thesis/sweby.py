import sys
import warnings
import numpy as np
import matplotlib as plt

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

DPI = 300
figformat = 'png'
fontsize = 16
plot_size = 10.0

torch.set_font_sizes(fontsize)

outputfile = 'sweby.' + figformat

#lines

r1 = np.arange(0.0, 1.0, 0.5)
y1 = 2.0 * r1

r2 = np.arange(0.0, 2.1, 0.5)
y2 = r2

r3 = np.arange(2.0, 4.1, 0.5)
y3 = np.empty(len(r3))
y3.fill(2.0)

r4 = np.arange(0.5, 4.1, 0.5)
y4 = np.empty(len(r4))
y4.fill(1.0)

fill1_r = [0, 1, 0.5]
fill1_y = [0, 1, 1]

fill2_r = [1, 4, 4, 2]
fill2_y = [1, 1, 2, 2]

albada_r = np.arange(0.0, 4.01, 0.01)
albada_y = (albada_r * albada_r + albada_r) / (albada_r * albada_r + 1.0)

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
asp_rat = 2.5/4.0
grid = plotter.axes1D((1,1), aspect_ratio=asp_rat)
grid[0].xaxis.set_ticks(np.arange(0, 4.1, 1.0))
grid[0].yaxis.set_ticks(np.arange(0, 3.1, 0.5))
grid[0].set_xlim([0.0, 4.0])
grid[0].set_ylim([0.0, 2.5])
grid[0].set_xlabel(plotter.format_label(torch.VarType('r', False)))
grid[0].set_ylabel(plotter.format_label(torch.VarType('\phi(r)', False)))

grid[0].xaxis.set_minor_locator(plt.ticker.FixedLocator([0.5, 1.5, 2.5, 3.5]))
grid[0].yaxis.set_minor_locator(plt.ticker.FixedLocator([0.25, 0.75, 1.25, 1.75, 2.25]))

### Plot.
kx = dict(linewidth=2.0, color='k')

grid[0].fill(fill1_r, fill1_y, color='0.5', alpha=0.5)
grid[0].fill(fill2_r, fill2_y, color='0.5', alpha=0.5)
grid[0].plot(albada_r, albada_y, linewidth=1.5, color='m')
grid[0].plot(r1, y1, **kx)
grid[0].plot(r2, y2, **kx)
grid[0].plot(r3, y3, **kx)
grid[0].plot(r4, y4, **kx)


###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile