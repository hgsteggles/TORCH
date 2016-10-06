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
plot_size = 3
fontsize = 10

torch.set_font_sizes(fontsize)

dir = "data/shadowing2/"
inputfile = []
inputfile.append(dir + "data2D_10.txt")
inputfile.append(dir + "data2D_25.txt")
inputfile.append(dir + "data2D_40.txt")
inputfile.append(dir + "data2D_55.txt")

outputfile = "multi-shadowing2" + '.' + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []
color_maps = []
for i in range(4):
	datacubes.append(torch.CFD_Data(inputfile[i], axial=False))
	vs_types.append(torch.VarType("nhii", isLog10=datacubes[i].appropriate_to_log("nhii")))
	vsminmax.append([2, 5])
	color_maps.append(hgspy.get_par_cmap(isReversed=False))

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, "linear", (2, 2), color_maps, tight=False, detail="all")

### Plotting.
plotter = torch.Plotter(datacubes[0].nx, datacubes[0].ny, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)

ts = 0.04
for i in range(len(datacubes)):
	timestring = str(datacubes[i].t) + " yrs"
	grid[i].text(1-ts, 1-ts, timestring, fontsize=int(1.5*fontsize), color='white', horizontalalignment='right', verticalalignment='top', transform = grid[i].transAxes)

for i in range(len(datacubes)):
	grid[i].xaxis.set_ticks(np.arange(0.03, 0.15 + 0.0001, 0.03))
	grid[i].yaxis.set_ticks(np.arange(0.03, 0.15 + 0.0001, 0.03))
	grid.cbar_axes[0].yaxis.set_ticks(np.arange(2.0, 5.0 + 0.0001, 0.5))

grid[2].xaxis.set_ticks(np.arange(0, 0.15 + 0.0001, 0.03))
grid[2].yaxis.set_ticks(np.arange(0, 0.15 + 0.0001, 0.03))

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch-multi in ' + outputfile

