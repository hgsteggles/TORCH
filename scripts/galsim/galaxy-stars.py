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

inputfile1 = "data/galsim/starpop-hm-a.txt"

data1 = np.genfromtxt(inputfile1, skip_header=1)

def plot_scatter(ax, dat, colormap):
    x = dat[:,12]
    y = dat[:,13]

    xy = np.vstack([x,y])
    z = scipy.stats.gaussian_kde(xy)(xy)

    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    z = np.log(z * len(x))

    im = ax.scatter(x, y, c=z, s=25, edgecolor='', cmap=colormap)

    ax.set_xlabel('$x\ \\left[\mathrm{kpc}\\right]$')
    ax.set_ylabel('$y\ \\left[\mathrm{kpc}\\right]$')

    return im

### Plotting
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotparams = torch.PlotParams(None, None, None, False, 'linear', (1, 1), None, tight=False, detail="all")
grid = plotter.getGrid(plotparams)
plotter.modifyGrid(grid, True)

ax = grid[0]

xline = [0, 20]
yline = [8.5, 8.5 - xline[1] / np.tan(np.radians(10.0))]
y2line = [8.5, 8.5 - xline[1] / np.tan(np.radians(65.0))]
ax.fill_between(xline, yline, y2line, facecolor='gray', alpha=0.5)

im = plot_scatter(ax, data1, hgspy.get_thermal_cmap())
ax.scatter([0], [8.5], c='k', marker='*', linewidths=[0.4], zorder=2)

ax.plot(xline, yline, color='gray')
ax.plot(xline, y2line, color='gray')

grid.cbar_axes[0].colorbar(im)
cbax = grid.cbar_axes[0]
cblax = cbax.axis[cbax.orientation]
#cbax.set_ylim([-0.1, 0.7])
cblax.label.set_text(plotter.format_label(torch.VarType('P', isLog10=True, units='kpc^{-2}')))

ax.set_xlim([-20, 20])
ax.set_ylim([-20, 20])

plotter.save_plot("galaxy-stars.png")