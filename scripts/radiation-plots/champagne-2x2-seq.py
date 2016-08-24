import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
import torch
load_src("hgspy", "../torchpack/hgspy.py")
import hgspy

DPI = 300
figformat = 'png'
plot_size = 10
fontsize = 32

torch.set_font_sizes(fontsize)

dir = "data/champagne/"
inputfile = []
inputfile.append(dir + "data2D_001.txt")
inputfile.append(dir + "data2D_025.txt")
inputfile.append(dir + "data2D_050.txt")
inputfile.append(dir + "data2D_100.txt")

outputfile = "multi-plot" + '.' + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []
color_maps = []
for i in range(4):
	datacubes.append(torch.CFD_Data(inputfile[i], axial=False))
	vs_types.append(torch.VarType("nhii", isLog10=datacubes[i].appropriate_to_log("nhii")))
	vsminmax.append([2, 4])
	color_maps.append(hgspy.get_par_cmap(isReversed=False))

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, "cubic", (2, 2), color_maps, tight=False, detail="all")
#plotparams.xminmax = [0, 0.5]
#plotparams.yminmax = [0.1, 0.6]

### Plotting.
plotter = torch.Plotter(1.0, 1.0, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)

ts = 0.04
for i in range(len(datacubes)):
	timestring = str(datacubes[i].t) + " yrs"
	grid[i].text(1-ts, 1-ts, timestring, fontsize=int(1.5*fontsize), color='white', horizontalalignment='right', verticalalignment='top', transform = grid[i].transAxes)

xminmax = [0, 0.5]
yminmax = [0.1, 0.6]
for i in range(len(datacubes)):
	grid[i].xaxis.set_ticks(np.arange(xminmax[0] + 0.1, xminmax[1] + 0.0001, 0.1))
	grid[i].yaxis.set_ticks(np.arange(yminmax[0] + 0.1, yminmax[1] + 0.0001, 0.1))
	grid[i].set_xlim(xminmax)
	grid[i].set_ylim(yminmax)
	grid.cbar_axes[0].yaxis.set_ticks(np.arange(2.0, 4.0 + 0.0001, 0.5))
grid[2].xaxis.set_ticks(np.arange(xminmax[0], xminmax[1] + 0.0001, 0.1))
grid[2].yaxis.set_ticks(np.arange(yminmax[0], yminmax[1] + 0.0001, 0.1))


###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch-multi in ' + outputfile

