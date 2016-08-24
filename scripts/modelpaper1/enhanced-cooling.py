import numpy as np
import warnings
import sys

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

DPI = 300
figformat = 'png'
plot_size = 5
fontsize = 16
outputfile = "enhanced-cooling.png"
torch.set_font_sizes(fontsize)

YR2S = 3.154e7
CM2PC = 3.24078e-19
PC2CM = 1.0 / CM2PC
snapshots = range(2, 52, 2)
data_id = 22
pad_data_id = "%02d" % (data_id,)
fileprefix = "data/model_paper1/set2/data_" + pad_data_id + "/data2D_"

### Data set up.
data = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)
mdot = data[:,8] # [g.s-1]
vinf = data[:,9] # [cm.s-1]

def getEnergyData(snapshot):
	pad_index = "%03d" % (snapshot,)
	datacube = torch.CFD_Data(fileprefix + pad_index + ".txt", axial=False)

	gamma = 5.0 / 3.0
	tem = datacube.get_var_raw('tem')
	thermal_energy = datacube.get_var_raw('pre') / (gamma - 1.0)
	kinetic_energy = datacube.get_var_raw('ke')
	vol = datacube.get_var_raw('vol-cylindrical') * PC2CM * PC2CM * PC2CM

	thermal_energy = thermal_energy * vol
	kinetic_energy = kinetic_energy * vol
	thermal_energy = np.delete(thermal_energy, np.where(tem < 1.0e5), 0)
	kinetic_energy = np.delete(kinetic_energy, np.where(tem < 1.0e5), 0)

	#return 100.0*np.sum(thermal_energy)
	#return 100.0*np.sum(kinetic_energy)
	return np.sum(thermal_energy + kinetic_energy)

def getEnergy(snapshot):
	j = (data_id - 1)%9
	return 2 * snapshot * 1000 * YR2S * 0.5 * mdot[j] * vinf[j] * vinf[j]

times = []
energies = []
energies_a = []
for i in range(len(snapshots)):
	times.append(snapshots[i] * 2)
	energies.append(getEnergyData(snapshots[i]))
	energies_a.append(getEnergy(snapshots[i]))

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
grid = plotter.axes1D((1,1), aspect_ratio=0.75)
grid[0].set_xlabel(plotter.format_label(torch.VarType('\\mathrm{Age}', units='kyr', isLog10=False)))
grid[0].set_ylabel(plotter.format_label(torch.VarType('E', units='erg', isLog10=True)))

### Plot.
grid[0].plot(times, energies, color='r', linewidth=1)
#grid[0].plot(times, energies_a, color='b', linewidth=1)
#grid[0].plot(times, np.log10(energies_a), label='Injected Energy', color='b', linestyle='--', linewidth=1)
#grid[0].plot(times, np.log10(energies), label='Total Energy', color='r', linewidth=1)
#grid[0].set_ylim([0, 1.0e44])

grid[0].legend(loc='center right')

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile