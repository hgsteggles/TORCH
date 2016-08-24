import math
import scipy.ndimage as ndimage
import matplotlib.ticker as ticker
import numpy as np
from astropy.io import fits

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData()

DPI = 300
figformat = 'png'
plot_size = 5
fontsize = 16
torch.set_font_sizes(fontsize=16)

hdu_list0 = fits.open(mp1_data.getRadioDirname(5, 2, 25, 45, 1.4) + "/intensity_beam_ff.fits")
hdu_list1 = fits.open(mp1_data.getRadioDirname(5, 2, 25, 45, 5) + "/intensity_beam_ff.fits")
data0 = hdu_list0[0].data
data1 = hdu_list1[0].data

sig = 5
FWHM = sig*hdu_list0[0].header['PIXAS']*2.335
data0 = ndimage.gaussian_filter(data0, sigma=(sig, sig), order=0)
data1 = ndimage.gaussian_filter(data1, sigma=(sig, sig), order=0)
#data0[data0 < 0.4] = 0.4
#data1[data1 < 0.4] = 0.4

sindex = (np.log10(data1 / 1000.0) - np.log10(data0 / 1000.0)) / ((np.log10(5.0) - np.log10(1.4)))

masked_array = np.ma.array (sindex, mask=np.logical_or(data0 < 0.4, data1 < 0.4))
cmap = hgspy.get_cdom_cmap(isReversed=True)
cmap.set_bad('w', 1.)

#smin = sindex[np.isfinite(sindex)].min()

#for x in np.nditer(sindex):
#	if not np.isfinite(x):
#		x = smin

nx = hdu_list0[0].header['NAXIS1']
ny = hdu_list0[0].header['NAXIS2']

### Plotting
plotter = torch.Plotter(nx, ny, plot_size, figformat, DPI)
plotparams = torch.PlotParams(None, None, None, False, 'linear', (1, 1), None, tight=False, detail="all")
grid = plotter.getGrid(plotparams)
plotter.modifyGrid(grid, True)

ax = grid[0]

dx = -hdu_list0[0].header['CDELT1'] * 60.0 * 60.0
dy = hdu_list0[0].header['CDELT2'] * 60.0 * 60.0

xmax = 0.5 * nx * dx
xmin = -0.5 * nx * dx
ymax = 0.5 * ny * dy
ymin = -0.5 * ny * dy

x = np.arange(xmin + (dx/2.0), xmax, dx)
y = np.arange(ymin + (dy/2.0), ymax, dy)
X, Y = np.meshgrid(x, y)
im = ax.imshow(masked_array, extent=[xmin,xmax,ymin,ymax], origin='lower', cmap=cmap)

formatter = ticker.ScalarFormatter(useOffset=True, useMathText=True)
formatter.set_powerlimits((0, 1))

grid.cbar_axes[0].colorbar(im)
cbax = grid.cbar_axes[0]
cblax = cbax.axis[cbax.orientation]
#cbax.set_ylim([-0.1, 0.7])
cblax.label.set_text(plotter.format_label(torch.VarType('\mathrm{Spectral\ Index}')))

x_range = 120.0

ax.set_xlim([-x_range/2.0, x_range/2.0])
ax.set_ylim([-x_range/2.0, x_range/2.0])
ax.xaxis.set_ticks(np.arange(-x_range/2.0, (x_range/2.0) + 0.0001, x_range/4.0))
ax.yaxis.set_ticks(np.arange(-(x_range/2.0), (x_range/2.0) + 0.0001, x_range/4.0))

ax.xaxis.get_major_ticks()[0].label1.set_visible(False)
ax.xaxis.get_major_ticks()[2].label1.set_visible(False)
ax.xaxis.get_major_ticks()[-1].label1.set_visible(False)

ax.yaxis.get_major_ticks()[0].label1.set_visible(False)
ax.yaxis.get_major_ticks()[2].label1.set_visible(False)
ax.yaxis.get_major_ticks()[-1].label1.set_visible(False)

plotter.save_plot("spectral-index-map.png")