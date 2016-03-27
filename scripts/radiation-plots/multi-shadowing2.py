import sys
import warnings
import numpy as np

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
dir = "data/shadowing2/"
fontsize = 16
var_type = "nhii"

torch.set_font_sizes(fontsize)

inputfile = []
inputfile.append(dir + "data2D_10.txt")
inputfile.append(dir + "data2D_25.txt")
inputfile.append(dir + "data2D_55.txt")

outputfile = 'multi-shadowing2' + '.' + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []
color_maps = []

cmap = hgspy.get_par_cmap(isReversed=False)

for i in range(3):
	datacubes.append(torch.CFD_Data(inputfile[i], axial=False))
	color_maps.append(cmap)
	vs_types.append(torch.VarType(var_type, datacubes[i].appropriate_to_log(var_type)))
	vsminmax.append([2.0, 5.0])

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, 'linear', (3, 1), color_maps, tight=False, detail="all")

### Plotting
plotter = torch.Plotter(datacubes[0].nx, datacubes[0].ny, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)

ts = 0.04
for i in range(len(grid)):
	timestring = str(datacubes[i].t) + " yrs"
	grid[i].text(1-ts, 1-ts, timestring, fontsize=int(1.5*fontsize), color='white', horizontalalignment='right', verticalalignment='top', transform = grid[i].transAxes)

for i in range(len(grid)):
	grid[i].xaxis.set_ticks(np.arange(0, 0.15 + 0.0001, 0.03))
	if i == 2:
		grid[i].yaxis.set_ticks(np.arange(0, 0.15 + 0.0001, 0.03))
	else:
		grid[i].yaxis.set_ticks(np.arange(0.03, 0.15 + 0.0001, 0.03))
	grid[i].xaxis.set_ticks_position('none')
	grid[i].yaxis.set_ticks_position('none')
	grid.cbar_axes[0].yaxis.set_ticks(np.arange(2.0, 5.0, 0.5))

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted ' + outputfile