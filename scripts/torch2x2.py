import argparse
import os
import sys
import warnings

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "./torchpack/torch.py")
load_src("hgspy", "./torchpack/hgspy.py")

import torch
import hgspy


DPI = 300
figformat = 'jpg'
plot_size = 10
fontsize = 20

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
color_maps = []
for i in range(4):
	datacubes.append(torch.CFD_Data(inputfile, axial=True))
	color_maps.append(hgspy.get_green_cmap())

vs_types = []
vs_types.append(torch.VarType("nh", datacubes[0].appropriate_to_log("nh")))
vs_types.append(torch.VarType("pre", datacubes[1].appropriate_to_log("pre")))
vs_types.append(torch.VarType("nhii", datacubes[2].appropriate_to_log("nhii")))
vs_types.append(torch.VarType("tem", datacubes[3].appropriate_to_log("tem")))

vsminmax = []
vsminmax.append([3.5, 4.8])
vsminmax.append([-10, -5])
vsminmax.append([2   , 4])
vsminmax.append([1, 5])

if args.minmax != None:
	dataminmax = torch.DataSetMinMax(args.minmax)
	vsminmax[0] = dataminmax.get('log10_nh')
	vsminmax[1] = dataminmax.get('log10_pre')
	vsminmax[2] = dataminmax.get('log10_tem')

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, "nearest", (2, 2), color_maps, tight=True, detail="some")

### Plotting
plotter = torch.Plotter(datacubes[0].nx, datacubes[0].ny, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)
for i in range(4):
	for axis in ['top','bottom','left','right']:
		grid[i].spines[axis].set_linewidth(0.5)
		grid.cbar_axes[i].spines[axis].set_linewidth(0.5)

#grid.cbar_axes[0].set_xlabel(r"$\mathregular{log_{10}(n_H}$ / $\mathregular{cm^{-3})}$")
#grid.cbar_axes[1].set_xlabel(r"$\mathregular{log_{10}(P}$ / $\mathregular{Ba)}$")
#grid.cbar_axes[2].set_xlabel(r"$\mathregular{log_{10}(n_{HII}}$ / $\mathregular{cm^{-3})}$")
#grid.cbar_axes[3].set_xlabel(r"$\mathregular{log_{10}(T}$ / $\mathregular{K})$")

plotter.add_quiver(5, datacubes[2], grid[2])

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

if args.s == False:
	print sys.argv[0] + ': plotted torch2x2 in ' + outputfile

