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

log = False
view_cbar = True
view_quiver = True
DPI = 300
figformat = 'jpg'
plot_size = 5

###	Parse arguements
parser = argparse.ArgumentParser(description='Plots 2D image of CFD data.')
parser.add_argument('inputfile', metavar='inputfile', type=str, help='Input file to produce image.')
parser.add_argument('var_type', metavar='var_type', type=str, choices=TorchCFD.var_types, help='Variable type.')
parser.add_argument('-v0', metavar='var_min', type=float, help='Minimum variable value.')
parser.add_argument('-v1', metavar='var_max', type=float, help='Maximum variable value.')
args = parser.parse_args()

inputfile = args.inputfile
outputfile = os.path.splitext(inputfile)[0] + '.' + figformat

###	Data set up.
torchData = TorchCFD(inputfile, axial=True)
var = torchData.get_var(args.var_type)
if log == True and torchData.appropriate_to_log(args.var_type):
	var = torchData.safe_log10(var)

var_min = args.v0 if args.v0 != None else var.min()
var_max = args.v1 if args.v1 != None else var.max()

### Plotting
plotter = TorchPlotter(torchData, plot_size, figformat, DPI)

###	Interpolation set up.
vari = torchData.interpolate(var, 'nearest')

###	Image.
ax = plotter.single(args.var_type, vari, var_min, var_max, hgspy.get_green_map())

###	Save figure.
plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch2x2 in ' + outputfile

