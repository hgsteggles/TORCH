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

params = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

massdat = params[:,0] # [Msun]
logLdat = params[:,2] # log[Q / phot.s-1]

logL_interp = interpolate.interp1d(massdat, logLdat)

print "check: " + str(logL_interp(9))

benstr = ""
if args.ben:
	benstr = "-ben"
if args.harry:
	benstr = "-harry"
outputfile = cornish_data.dirname + '/hist-lifetimes' + benstr + '.' + figformat

print benstr

### Data
simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + benstr + ".txt", skip_header=1)
cornish_survey = np.genfromtxt("data/cornish/cornish-uchiis.txt", skip_header=1, delimiter=',')

loglmin = 4.0
loglmax = 6.0
nbins = 10

bins = np.logspace(loglmin, loglmax, nbins + 1)
bincentres = np.sqrt(bins[:-1] * bins[1:])
logbins = np.linspace(loglmin, loglmax, nbins + 1)

print logbins

ages = []
for i in range(nbins):
	ages.append([])

for i in range(len(simulated_survey[:,0])):
	logL = logL_interp(simulated_survey[i, 1])

	if logL > logbins[0]:
		bin_index = 0

		for j in range(nbins):
			if logL <= logbins[j + 1]:
				bin_index = j
				break

		ages[bin_index].append(simulated_survey[i, 2])

for i in range(nbins):
	nremove = int(0.05 * len(ages[i]))
	#nremove = 0
	for j in range(nremove):
		del ages[i][0]
		del ages[i][-1]

phase_lifetimes = []

for i in range(nbins):
	if len(ages[i]) != 0:
		phase_lifetimes.append((np.max(ages[i]) - np.min(ages[i])))
		lifemax = np.max(ages[i])
	else:
		phase_lifetimes.append(0)
		lifemax = 0
	print str(len(ages[i])).rjust(3) + " " + str(int(phase_lifetimes[i])).rjust(6) + " " + str(int(lifemax)).rjust(6)

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)
plotter.ticklength *= 0.5

###	Axes.
asp_rat = 0.4
grid = plotter.axes1D((1,1), aspect_ratio=asp_rat)

grid[0].set_xlabel(plotter.format_label(torch.VarType('L', isLog10=True, units='L_\\odot')))
grid[0].set_ylabel(plotter.format_label(torch.VarType('\\mathrm{Phase Lifetimes}', units='yr')))
grid[0].set_xlim([loglmin, loglmax])
grid[0].set_xscale("log")
grid[0].set_yscale("log")
kx1 = dict(linewidth=1.5, label="Mottram et al. (2011)", color='b')
kx2 = dict(linewidth=1.5, label="Simulated", color='r', linestyle='--')

### Legend
handles, labels = grid[0].get_legend_handles_labels()
legend = grid[0].legend(handles, labels, loc=1)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile
