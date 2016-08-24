import math
import scipy.ndimage as ndimage
import matplotlib.ticker as ticker
import numpy as np
from astropy.io import fits
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
outputfile = "sed." + figformat
plot_size = 5
fontsize = 16
torch.set_font_sizes(fontsize)

fileprefix2 = "data/model_paper1/set2/data_33"

imagetype = "intensity_pixel_ff.fits"

def getFreqFlux(prefix, i):
	padi = "%02d" % ((i),)
	filesuffix = "/SED/radio_" + padi + "/" + imagetype
	hdu_list0 = fits.open(prefix + filesuffix)
	data0 = hdu_list0[0].data

	return hdu_list0[0].header['FREQ'] / 1.0e9, np.sum(data0) / 1000.0

freqs = []
fluxes = []
ufluxes = []
fluxes2 = []
ufluxes2 = []

for i in range(21):
	freq, flux = getFreqFlux("data/model_paper1/set2/data_24", i)
	freqs.append(freq)
	fluxes.append(flux)
	freq, flux = getFreqFlux("data/model_paper1/uniform", i)
	ufluxes.append(flux)
	freq, flux = getFreqFlux("data/model_paper1/set2/data_33", i)
	fluxes2.append(flux)
	freq, flux = getFreqFlux("data/model_paper1/uniform2", i)
	ufluxes2.append(flux)

#freqs = np.log10(np.array(freqs))
#fluxes = np.log10(np.array(fluxes))

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
grid = plotter.axes1D((1,1), aspect_ratio=0.67)
grid[0].set_xlabel(plotter.format_label(torch.VarType('\\nu', units='GHz')))
grid[0].set_ylabel(plotter.format_label(torch.VarType('S_{\\nu}', units='Jy')))

### Plot.
grid[0].loglog(freqs, fluxes, color='r', linewidth=1, label=r'$\alpha = 2$')
grid[0].loglog(freqs, ufluxes, color='b', linewidth=1, label=r'$\alpha = 0$')
#grid[0].loglog(freqs, fluxes2, color='g', linewidth=1, label=r'$\alpha = 2$')
#grid[0].loglog(freqs, ufluxes2, color='m', linewidth=1, label=r'$\alpha = 0$')

grid[0].set_xlim([0.3, 20])
grid[0].set_ylim([4, 100])

def myLogFormat(y,pos):
    decimalplaces = int(np.maximum(-np.log10(y),0))     # =0 for numbers >=1
    formatstring = '{{:.{:1d}f}}'.format(decimalplaces)
    return formatstring.format(y)

formatter = ticker.FuncFormatter(myLogFormat)
#formatter = ticker.ScalarFormatter(useOffset=True, useMathText=True)
#formatter.set_powerlimits((-2, 3))

grid[0].xaxis.set_major_formatter(formatter)
grid[0].yaxis.set_major_formatter(formatter)

### Legend
handles, labels = grid[0].get_legend_handles_labels()
legend = grid[0].legend(handles, labels, loc=4)

legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile


