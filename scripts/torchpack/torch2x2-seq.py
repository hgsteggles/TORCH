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

dir = "data/strickland2/"
inputfile = []
inputfile.append(dir + "data2D_15.txt")
inputfile.append(dir + "data2D_30.txt")
inputfile.append(dir + "data2D_45.txt")
inputfile.append(dir + "data2D_65.txt")

outputfile = "multi-plot" + '.' + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []
color_maps = []
for i in range(4):
	datacubes.append(torch.CFD_Data(inputfile[i], axial=False))
	vs_types.append(torch.VarType("nhii", datacubes[i].appropriate_to_log("nhii")))
	vsminmax.append([-1, 3])
	color_maps.append(hgspy.get_green_map())

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, "linear", (2, 2), color_maps, tight=False, detail="all")

### Plotting.
plotter = torch.Plotter(datacubes[0].nx, datacubes[0].ny, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)

### Extras.
for i in range(len(grid)):
	grid[i].xaxis.set_ticks(np.arange(0.5, 3.5, 0.5))
	grid[i].yaxis.set_ticks(np.arange(0.5, 3.5, 0.5))
	grid[i].xaxis.set_ticks_position('none')
	grid[i].yaxis.set_ticks_position('none')
grid[2].xaxis.set_ticks(np.arange(0, 3.5, 0.5))
grid[2].yaxis.set_ticks(np.arange(0, 3.5, 0.5))
grid.cbar_axes[0].yaxis.set_ticks(np.arange(-1.0, 3.5, 0.5))

ts = 0.04
for i in range(len(grid)):
	timestring = str(datacubes[i].t) + " yrs"
	grid[i].text(1-ts, 1-ts, timestring, fontsize=10, color='white', horizontalalignment='right', verticalalignment='top', transform = grid[i].transAxes)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch-multi in ' + outputfile

