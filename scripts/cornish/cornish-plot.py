import numpy as np
from astropy.io import fits
import matplotlib.ticker as ticker
import argparse
import scipy.ndimage as ndimage

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
import torch
load_src("hgspy", "../torchpack/hgspy.py")
import hgspy
load_src("fmt", "../torchpack/formatting.py")
import fmt
load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData()

DPI = 300
figformat = 'png'
plot_size = 5
fontsize = 10
torch.set_font_sizes(fontsize)

###	Parse arguements
parser = argparse.ArgumentParser(description='Plots CORNISH fits image.')
parser.add_argument('inputfile', metavar='inputfile', type=str, help='Input file to produce image.')
args = parser.parse_args()

inputfile = args.inputfile

hdu = fits.open(inputfile)[0]
nx = hdu.header['NAXIS1']
ny = hdu.header['NAXIS2']
dx = hdu.header['CDELT1']
dy = hdu.header['CDELT2']
nrefx = hdu.header['CRPIX1']
nrefy = hdu.header['CRPIX2']
vrefx = hdu.header['CRVAL1']
vrefy = hdu.header['CRVAL2']

imagexmin = vrefx + (0 - nrefx) * dx
imagexmax = vrefx + (nx - nrefx) * dx
imageymin = vrefy + (0 - nrefy) * dy
imageymax = vrefy + (ny - nrefy) * dy

xmin = imagexmin
xmax = imagexmax
ymin = imageymin
ymax = imageymax

### CASA?
CASA = False

if CASA:
	data = hdu.data[0][0]
	unitstr = 'mJy\,beam^{-1}'
	data = 1000.0 * data
else:
	data = hdu.data
	sig = 5
	FWHM = sig*hdu.header['PIXAS']*2.335
	data = ndimage.gaussian_filter(data, sigma=(sig, sig), order=0)
	#xmin = vrefx + (0 - nrefy) * dx
	#xmax = vrefx + (ny - nrefy) * dx
	unitstr = 'mJy\,pixel^{-1}'

### Plotting
plotter = torch.Plotter(nx, ny, plot_size, figformat, DPI)
plotparams = torch.PlotParams(None, None, None, False, 'linear', (1, 1), None, tight=False, detail="all")
grid = plotter.getGrid(plotparams)
plotter.modifyGrid(grid, True)

ax = grid[0]

im = ax.imshow(data, origin='lower', cmap='cubehelix', extent=[imagexmin,imagexmax,imageymin,imageymax])
ax.set_xlabel('RA')
ax.set_ylabel('Dec')
ax.set_xlim([xmin, xmax])
ax.set_ylim([ymin, ymax])
ax.set_xticks([xmin + 0.25*(xmax - xmin),
			  xmin + 0.5*(xmax - xmin),
			  xmin + 0.75*(xmax - xmin)])
ax.set_yticks([ymin + 0.25*(ymax - ymin),
			  ymin + 0.5*(ymax - ymin),
			  ymin + 0.75*(ymax - ymin)])

ax.xaxis.set_tick_params(color='w')
ax.yaxis.set_tick_params(color='w')
#ax.tick_params(axis='x', colors='white')
#ax.tick_params(axis='y', colors='white')

ra_formatter = ticker.FuncFormatter(fmt.fmt_ra)
dec_formatter = ticker.FuncFormatter(fmt.fmt_dec)
ax.xaxis.set_major_formatter(ra_formatter)
ax.yaxis.set_major_formatter(dec_formatter)

grid.cbar_axes[0].colorbar(im, format=ticker.FormatStrFormatter('%.2f'))
cbax = grid.cbar_axes[0]
cblax = cbax.axis[cbax.orientation]
cblax.label.set_text(plotter.format_label(torch.VarType('S_\\nu', isLog10=False, units=unitstr)))

plotter.save_plot("cornish-map.png")