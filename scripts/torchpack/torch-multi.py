import sys
import warnings
import matplotlib.pyplot as plt
import numpy as np
import torchpack.torch as torch
import torchpack.hgspy as hgspy



fontsize = 8
font = {'family':'STIXGeneral','style':'normal','variant':'normal','weight':'medium','size':(fontsize-2)}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':(fontsize+2)})
plt.rc('xtick',**{'labelsize':fontsize})
plt.rc('ytick',**{'labelsize':fontsize})
plt.rc('axes',**{'labelsize':(fontsize+2)})
#plt.rc('mathtext',**{'fontset':'custom','rm':'Bitstream Vera Sans','it':'Bitstream Vera Sans:italic','bf':'Bitstream Vera Sans:bold'})
plt.rc('mathtext',**{'fontset':'stix'})

DPI = 300
figformat = 'png'
plot_size = 5

inputfile = []
inputfile.append("data/shadowing2/data2D_10.txt")
inputfile.append("data/shadowing2/data2D_25.txt")
inputfile.append("data/shadowing2/data2D_40.txt")
inputfile.append("data/shadowing2/data2D_55.txt")

outputfile = "multi-plot" + '.' + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []
color_maps = []
for i in range(4):
	datacubes.append(torch.CFD_Data(inputfile[i], axial=False))
	vs_types.append(torch.VarType("nhii", datacubes[i].appropriate_to_log("nhii")))
	vsminmax.append([2, 5])
	color_maps.append(hgspy.get_green_map())

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, "nearest", (2, 2), color_maps, tight=False, detail="all")

### Plotting.
plotter = torch.Plotter(datacubes[0].nx, datacubes[0].ny, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)
for i in range(len(grid)):
	grid[i].xaxis.set_ticks(np.arange(0, 0.175, 0.025))
	grid[i].yaxis.set_ticks(np.arange(0, 0.175, 0.025))
	grid[i].xaxis.set_ticks_position('none')
	grid[i].yaxis.set_ticks_position('none')
grid[2].xaxis.set_ticks(np.arange(0, 0.15, 0.025))
grid[2].yaxis.set_ticks(np.arange(0, 0.15, 0.025))
grid.cbar_axes[0].yaxis.set_ticks(np.arange(2.2, 5.4, 0.4))

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch-multi in ' + outputfile

