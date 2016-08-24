import argparse
import os
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
plot_size = 10
fontsize = 32

torch.set_font_sizes(fontsize)

###	Parse arguements
parser = argparse.ArgumentParser(description='Plots 2D image of CFD data.')
parser.add_argument('-s', action='store_true', help='Run script silently.')
parser.add_argument('-minmax', metavar='minmaxfile', type=str, help='Minima and maxima of data.')
parser.add_argument('inputfile', metavar='inputfile', type=str, help='Input file to produce image.')
args = parser.parse_args()

inputfile = args.inputfile
outputfile = os.path.splitext(inputfile)[0] + '.' + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []

for i in range(4):
	datacubes.append(torch.CFD_Data(inputfile, axial=False))

color_maps = []
color_maps.append(hgspy.get_density_cmap())
color_maps.append(hgspy.get_water_cmap())
color_maps.append(hgspy.get_temperature_cmap())
color_maps.append(hgspy.get_par_cmap())


vs_types = []
vs_types.append(torch.VarType("nh", isLog10=datacubes[0].appropriate_to_log("nh")))
vs_types.append(torch.VarType("pre", isLog10=datacubes[1].appropriate_to_log("pre")))
vs_types.append(torch.VarType("tem", isLog10=datacubes[2].appropriate_to_log("tem")))
vs_types.append(torch.VarType("hii", isLog10=datacubes[3].appropriate_to_log("hii")))

vsminmax = []
#vsminmax.append([3.5, 4.8])
#vsminmax.append([-10, -5])
#vsminmax.append([2   , 4])
#vsminmax.append([1, 5])

vsminmax.append([None, None])
vsminmax.append([None, None])
vsminmax.append([None, None])
vsminmax.append([0, 1])

if args.minmax != None:
	dataminmax = torch.DataSetMinMax(args.minmax)
	vsminmax[0] = dataminmax.get('log10_nh')
	vsminmax[1] = dataminmax.get('log10_pre')
	vsminmax[2] = dataminmax.get('log10_tem')

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, "nearest", (2, 2), color_maps, tight=False, detail="all")

### Plotting
plotter = torch.Plotter(datacubes[0].nx * 0.5 / 0.6, datacubes[0].ny * 0.5 / 0.8, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)

plotter.add_quiver(grid[2], datacubes[2], 7, vmin=3.0e5, vmiddle=0.0)

xminmax = [0.0, 0.5]
yminmax = [0.1, 0.6]
for i in range(len(datacubes)):
	grid[i].xaxis.set_ticks(np.arange(xminmax[0], xminmax[1] + 0.0001, 0.1))
	grid[i].yaxis.set_ticks(np.arange(yminmax[0], yminmax[1] + 0.0001, 0.1))
	grid[i].set_xlim(xminmax)
	grid[i].set_ylim(yminmax)
	#grid.cbar_axes[0].yaxis.set_ticks(np.arange(2.0, 4.0 + 0.0001, 0.5))
#grid[2].xaxis.set_ticks(np.arange(xminmax[0], xminmax[1] + 0.0001, 0.1))
#grid[2].yaxis.set_ticks(np.arange(yminmax[0], yminmax[1] + 0.0001, 0.1))


###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

if args.s == False:
	print (sys.argv[0] + ': plotted torch2x2 in ' + outputfile)