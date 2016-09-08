import sys
import warnings
import numpy as np
from scipy import interpolate

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
parser.add_argument('-ben', action='store_true', help='Use simple ben stromgren survey.')
parser.add_argument('-harry', action='store_true', help='Use simple harry stromgren survey.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

load_src("cdat", "../torchpack/cornishdata.py")
import cdat
cornish_data = cdat.CornishData(args.iden)
star_data = cornish_data.star_data
nstars = len(star_data[:,0])

DPI = 300
figformat = 'png'
plot_size = 5.0
fontsize = 10
torch.set_font_sizes(fontsize)

benstr = ""
if args.ben:
	benstr = "-ben"
if args.harry:
	benstr = "-harry"
outputfile = cornish_data.dirname + '/hist-qlyc' + benstr + '.' + figformat

params = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

massdat = params[:,0] # [Msun]
logQdat = params[:,6] # log[Q / phot.s-1]

logQ_interp = interpolate.interp1d(massdat, logQdat)

### Data
simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + benstr + ".txt", skip_header=1)
cornish_qlyc = np.genfromtxt("data/cornish/Lyman_fluxes.txt", skip_header=1, dtype=float)[:,1]
cornish_qlyc = np.log10(cornish_qlyc)

qlyc = []
for i in range(len(simulated_survey[:,0])):
	qlyc.append(logQ_interp(simulated_survey[i,1]))

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotter.ticklength *= 0.5

###	Axes.
asp_rat = 0.4
grid = plotter.axes1D((1,1), aspect_ratio=asp_rat)

grid[0].set_xlabel(plotter.format_label(torch.VarType('\\log_{10}(Q_\mathrm{Lyc})', units='s^{-1}')))
grid[0].set_ylabel(plotter.format_label(torch.VarType('N')))
grid[0].set_xlim([43, 50])
kx1 = dict(linewidth=1.5, label="CORNISH", color='b')
kx2 = dict(linewidth=1.5, label="Simulated", color='r', linestyle='--')

bins = np.linspace(43.0, 50.0, 20)
bincentres = np.sqrt(bins[:-1] * bins[1:])

plotter.histstep(grid[0], cornish_qlyc, bins, errorcentres=bincentres, **kx1)
plotter.histstep(grid[0], qlyc, bins, **kx2)

### Legend
handles, labels = grid[0].get_legend_handles_labels()
legend = grid[0].legend(handles, labels, loc=1)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile
