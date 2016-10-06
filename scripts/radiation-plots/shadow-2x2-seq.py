import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np
import math

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

dir = "data/shadowing-square/"
inputfile = []
inputfile.append(dir + "data2D_002.txt")
inputfile.append(dir + "data2D_011.txt")
inputfile.append(dir + "data2D_040.txt")
inputfile.append(dir + "data2D_100.txt")

outputfile = "multi-shadow" + '.' + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []
color_maps = []
for i in range(4):
	datacubes.append(torch.CFD_Data(inputfile[i], axial=True))
	vs_types.append(torch.VarType("hii", isLog10=datacubes[i].appropriate_to_log("hii")))
	vsminmax.append([0, 1])
	color_maps.append(hgspy.get_par_cmap(isReversed=False))

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, "linear", (2, 2), color_maps, tight=False, detail="all")

### Plotting.
plotter = torch.Plotter(datacubes[0].nx, datacubes[0].ny, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)

PC2CM = 3.09e18
alphaB = 2.59e-13
srcS = 1.0e48
nh = datacubes[0].get_var_raw('nh')[0]
invtrec = nh*alphaB
trec = 1.0/invtrec

clumpx = 3*np.array([-0.1, -0.1, 0.1, 0.1, -0.1])
clumpy = 3*np.array([0.6, 0.8, 0.8, 0.6, 0.6])

RS_inf = math.pow((3.0*srcS)/(4.0*math.pi*nh*nh*alphaB), 1.0/3.0)/PC2CM

print RS_inf

RS = []
RS.append(RS_inf*math.pow(1.0 - math.exp(-datacubes[0].tsecs/trec), 1.0/3.0))
RS.append(RS_inf*math.pow(1.0 - math.exp(-datacubes[1].tsecs/trec), 1.0/3.0))
RS.append(RS_inf*math.pow(1.0 - math.exp(-datacubes[2].tsecs/trec), 1.0/3.0))
RS.append(RS_inf*math.pow(1.0 - math.exp(-datacubes[3].tsecs/trec), 1.0/3.0))

theta = np.linspace(0.0, 2.01*math.pi, 300, endpoint=True)
circx = []
circy = []
for i in range(4):
    circx.append(RS[i]*np.cos(theta))
    circy.append(1.5 + RS[i]*np.sin(theta))
circx_inf =  RS_inf*np.cos(theta)
circy_inf = 1.5 + RS_inf*np.sin(theta)

kwargs = {'linewidth':1}

for i in range(len(datacubes)):
	grid[i].plot(clumpx, clumpy, c='r', **kwargs)
	grid[i].plot(circx[i], circy[i], c='m', **kwargs)
	grid[i].plot(circx_inf, circy_inf, c='w', **kwargs)

for i in range(len(datacubes)):
	grid[i].xaxis.set_ticks(np.arange(-1.0, 1.5 + 0.0001, 0.5))
	grid[i].yaxis.set_ticks(np.arange(0.5, 3.0 + 0.0001, 0.5))
grid[2].xaxis.set_ticks(np.arange(-1.5, 1.5 + 0.0001, 0.5))
grid[2].yaxis.set_ticks(np.arange(0.0, 3.0 + 0.0001, 0.5))

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch-multi in ' + outputfile

