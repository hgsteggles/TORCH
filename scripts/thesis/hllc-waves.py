import sys
import warnings
import numpy as np
import matplotlib as plt


def load_src(name, fpath):
	import os, imp
	return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))


load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

DPI = 300
figformat = 'png'
fontsize = 16
plot_size = 10.0

torch.set_font_sizes(fontsize)

outputfile = 'hllc-waves.' + figformat


def arrowed_spines(ax=None, arrow_length=20, labels=('', ''), arrowprops=None):
	xlabel, ylabel = labels
	if ax is None:
		ax = plt.gca()
	if arrowprops is None:
		arrowprops = dict(arrowstyle='<|-', facecolor='black')

	for i, spine in enumerate(['left', 'bottom']):
		# Set up the annotation parameters
		t = ax.spines[spine].get_transform()
		xy, xycoords = [1, 0], ('axes fraction', t)
		xytext, textcoords = [arrow_length, 0], ('offset points', t)
		ha, va = 'left', 'bottom'

		# If axis is reversed, draw the arrow the other way
		top, bottom = ax.spines[spine].axis.get_view_interval()
		if top > bottom:
			xy[0] = 0
			xytext[0] *= -1
			ha, va = 'right', 'top'

		if spine is 'bottom':
			xarrow = ax.annotate(xlabel, xy, xycoords=xycoords, xytext=xytext,
								 textcoords=textcoords, ha=ha, va='center',
								 arrowprops=arrowprops)
		else:
			yarrow = ax.annotate(ylabel, xy[::-1], xycoords=xycoords[::-1],
								 xytext=xytext[::-1], textcoords=textcoords[::-1],
								 ha='center', va=va, arrowprops=arrowprops)
	return xarrow, yarrow  # lines

r1 = np.arange(0.0, 1.1, 0.5)
y1 = r1

r2 = np.arange(-0.5, 0.1, 0.5)
y2 = - 2.0 * r2

r3 = np.arange(0.0, 0.26, 0.25)
y3 = 4.0 * r3

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
asp_rat = 2.5 / 4.0
grid = plotter.axes1D((1, 1), aspect_ratio=asp_rat)
grid[0].set_xlim([-0.7, 1.2])
grid[0].set_ylim([0.0, 1.2])
grid[0].xaxis.set_ticks(np.array([0]))
grid[0].yaxis.set_ticks([])
grid[0].spines['left'].set_position('zero')
grid[0].spines['right'].set_color('none')
grid[0].spines['bottom'].set_position('zero')
grid[0].spines['top'].set_color('none')
grid[0].xaxis.set_ticks_position('bottom')
grid[0].yaxis.set_ticks_position('left')

arrprops = dict(arrowstyle='<|-, head_width=0.4, head_length=0.8',facecolor='black',linewidth=plotter.linewidth)
arrowed_spines(grid[0], arrowprops=arrprops)

# grid[0].set_xlabel(plotter.format_label(torch.VarType('x', False)))
# grid[0].set_ylabel(plotter.format_label(torch.VarType('t', False)))

# sets axes labels on both ends
grid[0].annotate('$x$', xy=(0.98, -0.01), ha='left', va='top',
				 xycoords='axes fraction', fontsize=20)
grid[0].annotate('$t$', xy=(0.37, 0.98), xytext=(-15, 2), ha='left', va='top',
				 xycoords='axes fraction', textcoords='offset points', fontsize=20)

grid[0].annotate('$\\vec{U}_\\mathrm{L}$', xy=(0.10, 0.5), ha='left', va='top',
				 xycoords='axes fraction', fontsize=24)
grid[0].annotate('$\\vec{U}_{\\star \\mathrm{L}}$', xy=(0.28, 0.5), ha='left', va='top',
				 xycoords='axes fraction', fontsize=24)
grid[0].annotate('$\\vec{U}_{\\star \\mathrm{R}}$', xy=(0.52, 0.5), ha='left', va='top',
				 xycoords='axes fraction', fontsize=24)
grid[0].annotate('$\\vec{U}_\\mathrm{R}$', xy=(0.74, 0.5), ha='left', va='top',
				 xycoords='axes fraction', fontsize=24)
grid[0].annotate('$S_\\mathrm{L}$', xy=(0.07, 0.91), ha='left', va='top',
				 xycoords='axes fraction', fontsize=24)
grid[0].annotate('$S_\\mathrm{R}$', xy=(0.94, 0.90), ha='right', va='top',
				 xycoords='axes fraction', fontsize=24)
grid[0].annotate('$S_{\star}$', xy=(0.53, 0.91), ha='right', va='top',
				 xycoords='axes fraction', fontsize=24)

### Plot.
kx = dict(linewidth=2.0, color='r')
kx2 = dict(linewidth=2.0, linestyle='--', color='r')

grid[0].plot(r1, y1, **kx)
grid[0].plot(r2, y2, **kx)
grid[0].plot(r3, y3, **kx2)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted in ' + outputfile
