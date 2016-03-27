import argparse
import os
import sys
import matplotlib.pyplot as plt
import torchpack.torch as torch
import torchpack.hgspy as hgspy



font = {'family':'STIXGeneral','style':'normal','variant':'normal','weight':'medium','size':18}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':18})
plt.rc('xtick',**{'labelsize':14})
plt.rc('ytick',**{'labelsize':14})
plt.rc('axes',**{'labelsize':18})
#plt.rc('mathtext',**{'fontset':'custom','rm':'Bitstream Vera Sans','it':'Bitstream Vera Sans:italic','bf':'Bitstream Vera Sans:bold'})
plt.rc('mathtext',**{'fontset':'stix'})

log = False
view_cbar = True
view_quiver = True
DPI = 300
figformat = 'png'
plot_size = 5

###	Parse arguements
parser = argparse.ArgumentParser(description='Plots 2D image of CFD data.')
parser.add_argument('inputfile', metavar='inputfile', type=str, help='Input file to produce image.')
parser.add_argument('var_type', metavar='var_type', type=str, choices=torch.CFD_Data.var_typenames, help='Variable type.')
parser.add_argument('-v0', metavar='var_min', type=float, help='Minimum variable value.')
parser.add_argument('-v1', metavar='var_max', type=float, help='Maximum variable value.')
args = parser.parse_args()

inputfile = args.inputfile
outputfile = os.path.splitext(inputfile)[0] + '.' + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []
color_maps = []

datacubes.append(torch.CFD_Data(inputfile, axial=False))
vs_types.append(torch.VarType(args.var_type, datacubes[0].appropriate_to_log(args.var_type)))
vsminmax.append([args.v0, args.v1])
color_maps.append(hgspy.get_green_map())

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, False, (1, 1), color_maps, tight=True, detail="all")

### Plotting
plotter = torch.Plotter(datacubes[0].nx, datacubes[0].ny, plot_size, figformat, DPI)

###	Image.
grid = plotter.multi(plotparams)

###	Save figure.
plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch1 in ' + outputfile

