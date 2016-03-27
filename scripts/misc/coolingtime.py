import argparse
import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from packages import hgspy
from torch_plot import TorchCFD, TorchCool, TorchPlotter



font = {'family':'serif','size':6}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':6})
plt.rc('xtick',**{'labelsize':6})
plt.rc('ytick',**{'labelsize':6})

fformat = 'jpg'
log = True

###	Parse arguements
parser = argparse.ArgumentParser(description='Plots 2D image of CFD cooling data.')
parser.add_argument('cfdfile', metavar='cfdfile', type=str, help='CFD file to produce image.')
parser.add_argument('coolfile', metavar='coolfile', type=str, help='Cooling file to produce image.')
args = parser.parse_args()

outputfile = os.path.splitext(args.cfdfile)[0] + '.' + fformat

###	Data set up.
cfd = TorchCFD(args.cfdfile, axial=True)
cool = TorchCool(args.coolfile, axial=True)

nh = cfd.get_var('nh')
tem = cfd.get_var('tem')
crate = cool.get_var('cool')
L = np.absolute(cool.get_var('lcool')/(nh*nh))
E = 1.5*nh*1.3806488e-16*tem
with np.errstate(divide='ignore', over='ignore'):
	dt = E / np.absolute(crate) / 3.15569e7
dt[crate == 0] = 1e10
dt[~np.isfinite(dt)] = 1e10
L[dt > 20] = L[dt < 20].min()

vs = [nh, tem, L, dt]
for i in range(3):
	vs[i] = cfd.safe_log10(vs[i])

vsminmax = [[-0.87, 6.89],
			[None, None],
			[None, None],
			[0, 20]]

### Plotting
plotter = TorchPlotter(cfd, plot_size=5, figformat=fformat, dpi=300)

###	Interpolation set up.
vsi = []
for i in range(4):
	vsi.append(cfd.interpolate(vs[i], 'nearest'))
	vsminmax[i][0] = (vs[i].min() if vsminmax[i][0] == None else vsminmax[i][0])
	vsminmax[i][1] = (vsi[i].max() if vsminmax[i][1] == None else vsminmax[i][1])

###	Image.
vs_types = ['nH', 'energy', 'lambda', 'dtcool']
grid = plotter.twobytwo(vs_types, vsi, vsminmax, hgspy.get_green_map(), detail="all")
for i in range(4):
	for axis in ['top','bottom','left','right']:
		grid[i].spines[axis].set_linewidth(0.5)
		grid.cbar_axes[i].spines[axis].set_linewidth(0.5)
grid.cbar_axes[0].set_xlabel(r"$\mathregular{log_{10}(n_H}$ / $\mathregular{cm^{-3})}$")
grid.cbar_axes[1].set_xlabel(r"$\mathregular{log_{10}(E}$ / $\mathregular{ergs.cm^{-3})}$")
grid.cbar_axes[2].set_xlabel(r"$\mathregular{log_{10}(\Lambda}$ / $\mathregular{ergs.s^{-1}.cm^{3})}$")
grid.cbar_axes[3].set_xlabel(r"$\mathregular{dt}$ / $\mathregular{yrs}$")

###	Save figure.
plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted torch2x2 in ' + outputfile

