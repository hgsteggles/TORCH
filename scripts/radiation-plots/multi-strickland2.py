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
dir = "data/strickland2/"
fontsize = 16
var_type = "nh"

torch.set_font_sizes(fontsize)

inputfile = []
inputfile.append(dir + "data2D_15.txt")
inputfile.append(dir + "data2D_45.txt")
inputfile.append(dir + "data2D_65.txt")

outputfile = 'multi-strickland2' + '.' + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []
color_maps = []

cmap = hgspy.get_density_cmap(isReversed=False)

for i in range(3):
	datacubes.append(torch.CFD_Data(inputfile[i], axial=False))
	color_maps.append(cmap)
	vs_types.append(torch.VarType(var_type, isLog10=datacubes[i].appropriate_to_log(var_type)))
	vsminmax.append([-1.0, 3.0])

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, 'linear', (3, 1), color_maps, tight=False, detail="all")

### Plotting
plotter = torch.Plotter(datacubes[0].nx, datacubes[0].ny, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)

ts = 0.04
for i in range(len(datacubes)):
	timestring = str(datacubes[i].t) + " yrs"
	grid[i].text(1-ts, 1-ts, timestring, fontsize=int(1.5*fontsize), color='white', horizontalalignment='right', verticalalignment='top', transform = grid[i].transAxes)

for i in range(len(datacubes)):
	grid[i].xaxis.set_ticks(np.arange(0, 3.0 + 0.0001, 0.5))
	grid[i].yaxis.set_ticks(np.arange(0, 3.0 + 0.0001, 0.5))
	grid[i].xaxis.set_ticks_position('none')
	grid[i].yaxis.set_ticks_position('none')
	grid.cbar_axes[0].yaxis.set_ticks(np.arange(-1.0, 3.0 + 0.00001, 0.5))

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted ' + outputfile