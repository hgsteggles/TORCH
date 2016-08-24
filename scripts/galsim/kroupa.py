import math
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
outputfile_qlyc = "kroupa.png"

torch.set_font_sizes(fontsize)

### Kroupa
p = []
p.append(0.3)
p.append(1.3)
p.append(2.3)

mlim = [0.01, 0.08, 0.5, 120]

A = [1, 1, 1]
A0 = 0

for i in range(3):
	for j in range(i):
		A[i] *= math.pow(mlim[j+1], p[j+1] - p[j])
	A0 += A[i] * (math.pow(mlim[i+1], 1 - p[i]) - math.pow(mlim[i], 1 - p[i])) / (1 - p[i])

for i in range(3):
	A[i] /= A0

def imf(m):
	for i in range(1, len(mlim)):
		if m <= mlim[i]:
			return A[i-1] * math.pow(m, -p[i-1])

### Data
N = 101
m = []
e = []

for i in range(len(mlim)-1):
	for x in range(0, N):
		m.append(mlim[i] + (mlim[i+1] - mlim[i]) * float(x) / float(N - 1))
		e.append(imf(m[-1]))

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
grid = plotter.axes1D((1,1), aspect_ratio=0.75)
grid[0].set_xlabel(plotter.format_label(torch.VarType('m', units='M_{\odot}', isLog10=False)))
grid[0].set_ylabel(plotter.format_label(torch.VarType('\\xi(m) \\Delta m', units='', isLog10=False)))

### Plot.
grid[0].loglog(m, e, color='m')

grid[0].set_xlim([0.01, 1000])
grid[0].set_ylim([min(e), 100])

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile_qlyc)

print sys.argv[0] + ': plotted in ' + outputfile_qlyc