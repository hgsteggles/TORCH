import linecache
import numpy as np
import scipy
import scipy.stats
import matplotlib.pyplot as plt

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
torch.set_font_sizes(fontsize=16)

inputfile1 = "data/galsim/starpop-hm-w-a.txt"

data1 = np.genfromtxt(inputfile1, skip_header=1)

def plot_scatter(ax, dat, colormap):
    x = dat[:,12]
    y = dat[:,13]

    xy = np.vstack([x,y])
    z = scipy.stats.gaussian_kde(xy)(xy)

    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    z = np.log(z)

    ax.scatter(x, y, c=z, s=25, edgecolor='', cmap=colormap)

    ax.set_xlabel('$x\ \\left[\mathrm{kpc}\\right]$')
    ax.set_ylabel('$y\ \\left[\mathrm{kpc}\\right]$')

### Plotting
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotparams = torch.PlotParams(None, None, None, False, 'linear', (1, 1), None, tight=False, detail="all")
grid = plotter.getGrid(plotparams)
plotter.modifyGrid(grid, False)

ax = grid[0]

plot_scatter(ax, data1, hgspy.get_density_cmap())

ax.set_xlim([-20, 20])
ax.set_ylim([-20, 20])

plotter.save_plot("galaxy-stars.png")