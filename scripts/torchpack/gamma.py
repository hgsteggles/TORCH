import argparse
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import math

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

import torch
import hgspy

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
fontsize = 16

torch.set_font_sizes(fontsize)

KB = 1.38e-16
RGAS = 8.3144598e7

###	Parse arguements
parser = argparse.ArgumentParser(description='Plots 2D image of CFD data.')
parser.add_argument('inputfile', metavar='inputfile', type=str, help='Input cooling file to produce image.')
parser.add_argument('cfdfile', metavar='cfdfile', type=str, help='Input cfd file to produce image.')
parser.add_argument('var_type', metavar='var_type', type=str, choices=torch.Cool_Data.var_typenames, help='Variable type.')
parser.add_argument('-v0', metavar='var_min', type=float, help='Minimum variable value.')
parser.add_argument('-v1', metavar='var_max', type=float, help='Maximum variable value.')
args = parser.parse_args()

inputfile = args.inputfile
cfdfile = args.cfdfile
outputfile = "gamma1." + figformat

###	Data set up.
datacubes = []
vs_types = []
vsminmax = []
color_maps = []

cfdcube = torch.CFD_Data(cfdfile, axial=True)
datacubes.append(torch.Cool_Data(inputfile, axial=True))

vs_types.append(torch.VarType(args.var_type, units="", isLog10=True))
vsminmax.append([args.v0, args.v1])
color_maps.append(hgspy.get_water_cmap())

plotparams = torch.PlotParams(datacubes, vs_types, vsminmax, True, "nearest", (1, 1), color_maps, tight=True, detail="all")

### Plotting
plotter = torch.Plotter(datacubes[0].nx, datacubes[0].ny, plot_size, figformat, DPI)

###	Image.
vs = datacubes[0].get_var_raw('cool')
vs *= np.sqrt(cfdcube.get_var_raw('tem')) / (cfdcube.get_var_raw('nh') * cfdcube.get_var_raw('nh'))

g = 5.0 / 3.0
muh = 1.0e-24
x = 1.0
C1 = (1.0 / vs) * np.sqrt(8.0 * (muh * (g - 1))**5 / (KB * x * (g + 1)**8))
vinf = 311000000.0
mdot = 9.79e18
L = 0.5 * mdot * vinf * vinf
rho_0 = cfdcube.get_var_raw('den')
vcrit = np.power(rho_0 * L / (2.0 * math.pi * (C1 * C1)), 1.0 / 11.0)

grid = plotter.multi(plotparams, vcrit / 1.0e5)

###	Save figure.
plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch1 in ' + outputfile

