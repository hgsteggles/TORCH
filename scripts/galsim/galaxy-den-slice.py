import linecache
import numpy as np
import scipy
import matplotlib.pyplot as plt
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
plot_size = 10
fontsize = 32
torch.set_font_sizes(fontsize)

inputfile = "data/galsim/density.txt"

### Read data.
data = np.genfromtxt(inputfile, skip_header=1)

nx = int(math.sqrt(len(data[:,0])))
ny = int(len(data[:,0]) / nx)

x = []
x.append(data[:,0])
x.append(data[:,1])

dx = data[1,1] - data[0,1]

xi = []
xi.append(np.linspace(x[0].min(), x[0].max(), nx))
xi.append(np.linspace(x[1].min(), x[1].max(), ny))
xi[0], xi[1] = np.meshgrid(xi[0], xi[1])

den = data[:,4]

deni = scipy.interpolate.griddata((x[0], x[1]), den,
						   (xi[0], xi[1]),
						   method="linear")

### Plotting
plotter = torch.Plotter(nx, ny, plot_size, figformat, DPI)
plotparams = torch.PlotParams(None, None, None, False, 'linear', (1, 1), None, tight=False, detail="all")
grid = plotter.getGrid(plotparams)
plotter.modifyGrid(grid, True)

ax = grid[0]
cbax = grid.cbar_axes[0]

var = deni

im = ax.imshow(var, vmin=var.min(), vmax=var.max(), origin='lower',
			   extent=[x[0].min(), x[0].max(), x[1].min(), x[1].max()],
			   interpolation="bicubic", cmap=hgspy.get_par_cmap())

ax.set_xlabel('$x\ \\left[\mathrm{kpc}\\right]$')
ax.set_ylabel('$y\ \\left[\mathrm{kpc}\\right]$')

cbax.colorbar(im)
cblax = cbax.axis[cbax.orientation]
#cbax.set_ylim([-0.1, 0.7])
cblax.label.set_text(plotter.format_label(torch.VarType('n_\mathrm{e}', units='cm^{-3}')))

plotter.save_plot("galaxy-density.png")

### Plotting
#plotter = torch.Plotter(nx, ny, plot_size, figformat, DPI)
#plotparams = torch.PlotParams(None, None, None, False, 'linear', (1, 1), None, tight=False, detail="all")
#grid = plotter.getGrid(plotparams)
#plotter.modifyGrid(grid, True)

#ax = grid[0]
