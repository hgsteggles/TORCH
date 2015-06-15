from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import hgspy
from torch_plot import TorchCFD, TorchInterpolatedGrid, TorchPlotter
import argparse
import os
import sys

font = {'family':'serif','size':6}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':6})
plt.rc('xtick',**{'labelsize':6})
plt.rc('ytick',**{'labelsize':6})

log = True
view_cbar = True
view_quiver = True
DPI = 300
figformat = 'jpg'
plot_size = 5

###	Parse arguements
parser = argparse.ArgumentParser(description='Plots 2D image of CFD data.')
parser.add_argument('inputfile', metavar='inputfile', type=str, help='Input file to produce image.')
args = parser.parse_args()

inputfile = args.inputfile
outputfile = os.path.splitext(inputfile)[0] + '.' + figformat

###	Data set up.
torchData = TorchCFD(inputfile, axial=True)
vs_types = ["nh", "pre", "nhii", "tem"]
vs = []
for i in range(4):
	vs.append(torchData.get_var(vs_types[i]))
	if log == True and torchData.appropriate_to_log(vs_types[i]):
		vs[i] = torchData.safe_log10(vs[i])

vsminmax = [[-3.08, 6.49],
			[-12.59, -7.04],
			[2   , 4   ],
			[0.49, 8.15]]

### Plotting
plotter = TorchPlotter(torchData, plot_size, figformat, DPI)

###	Interpolation set up.
vsi = []
for i in range(4):
	vsi.append(torchData.interpolate(vs[i], 'nearest'))
	vsminmax[i][0] = (vs[i].min() if vsminmax[i][0] == None else vsminmax[i][0])
	vsminmax[i][1] = (vsi[i].max() if vsminmax[i][1] == None else vsminmax[i][1])

###	Image.
grid = plotter.twobytwo(vs_types, vsi, vsminmax, hgspy.get_green_map(), detail="some")
for i in range(4):
	for axis in ['top','bottom','left','right']:
		grid[i].spines[axis].set_linewidth(0.5)
		grid.cbar_axes[i].spines[axis].set_linewidth(0.5)
grid.cbar_axes[0].set_xlabel(r"$\mathregular{log_{10}(n_H}$ / $\mathregular{cm^{-3})}$")
grid.cbar_axes[1].set_xlabel(r"$\mathregular{log_{10}(P}$ / $\mathregular{Ba)}$")
grid.cbar_axes[2].set_xlabel(r"$\mathregular{log_{10}(n_{HII}}$ / $\mathregular{cm^{-3})}$")
grid.cbar_axes[3].set_xlabel(r"$\mathregular{log_{10}(T}$ / $\mathregular{K})$")

###	Save figure.
plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch2x2 in ' + outputfile

