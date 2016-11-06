import warnings

import numpy as np
import sys
import math

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
fontsize = 12
torch.set_font_sizes(fontsize)

outputfile = 'dtype-1d-late-2000' + '.' + figformat

inputfile = "data/D-type/late-1D-hr-t-r/if.csv"

def Euler(ta, tb, n, fa, f):
	R = np.zeros(n + 1)
	t = np.zeros(n + 1)
	R[0] = fa
	t[0] = ta
	h = (tb - ta) / float(n)
	for i in range(n):
		R[i + 1] = (R[i] + h * f(t[i], R[i]))
		t[i + 1] = (t[i] + h)
	return R, t

###	Data set up.
S2YR = 1.0 / (365.0 * 24.0 * 60.0 * 60.0)
PC2CM = 3.09e18
RGAS = 8.3144598e7

data = np.genfromtxt(inputfile, delimiter=',')

t = data[:,0] / S2YR #[yrs]
R = data[:,1] * PC2CM #[cm]

nh = 3112.3
alphaB = 2.59e-13
Q = 1.0e49
RS = np.power(3.0 * Q / (4.0 * math.pi * nh * nh * alphaB), 1.0 / 3.0)

def soundSpeed(T, mu):
	return math.sqrt(RGAS * T / mu)

T_o = 1.0e3
T_i = 1.0e4
mu_o = 1.0
mu_i = 0.5
co = soundSpeed(T_o, mu_o)
#ci = soundSpeed(T_i, mu_i)
ci = 12.85e5
cico = ci / co
ts = RS / ci

def spitzer_exact(t, y):
	return ci * (RS / y)**(3.0 / 4.0) - ci * (mu_i * T_o / (mu_o * T_i)) * (RS / y)**(- 3.0 / 4.0)

def bisbas_exact(t, y):
	inner = (4.0 / 3.0) * (RS / y)**(3.0 / 2.0) - (mu_i * T_o / (2.0 * mu_o * T_i))
	inner = max(inner, 0)

	return ci * math.sqrt(inner)

def raga_exact(t, y):
	inner = (4.0 / 3.0) * (RS / y)**(3.0 / 2.0) - (mu_i * T_o / (mu_o * T_i))
	inner = max(inner, 0)

	return ci * math.sqrt(inner)

res = 1000

R_spit, t_spit = Euler(0, t[-1], res * t.shape[0] - 1, RS, spitzer_exact)
R_raga, t_raga = Euler(0, t[-1], res * t.shape[0] - 1, RS, raga_exact)
R_bisbas, t_bisbas = Euler(0, t[-1], res * t.shape[0] - 1, RS, bisbas_exact)
A = 0.733
f_SB = 1.0 - A * np.exp(-t_bisbas * S2YR / 1.0e6)
R_SB = R_bisbas + f_SB * (R_spit - R_bisbas)

error_SB = np.abs((R_SB[::res] - R) / R_SB[::res])

stag = (cico)**(4.0 / 3.0) * RS
R_stag = np.array([stag, stag])
R_stag2 = (4.0 / 3.0)**(2.0 / 3.0) * R_stag
t_stag = np.array([0, t[-1]])

### Plotting.
plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

###	Axes.
asp_rat = 7.0 / 16.0
grid = plotter.axes1D((2,1), aspect_ratio=asp_rat)
grid[0].yaxis.set_ticks(np.arange(0.005, 0.035, 0.005))
grid[1].yaxis.set_ticks(np.arange(0.0, 6, 1.0))
grid[0].set_xlim([0.0, 3.0])
grid[1].set_ylim([0, 5])
grid[1].set_xlabel(plotter.format_label(torch.VarType('t', units='Myr')))
grid[0].set_ylabel(plotter.format_label(torch.VarType('\\mathrm{Relative\ Error}')))
grid[1].set_ylabel(plotter.format_label(torch.VarType('R_\mathrm{IF}', units='pc')))

### Plot.
kx = dict(linewidth=1.2)
kx2 = dict(linewidth=1.0)

grid[1].plot(t * S2YR / 1.0e6, R / PC2CM, color='black', linestyle='-', label='Simulation', **kx)
grid[1].plot(t_spit * S2YR / 1.0e6, R_spit / PC2CM, color='blue', linestyle='-', label='Raga-I', **kx)
grid[1].plot(t_raga * S2YR / 1.0e6, R_raga / PC2CM, color='red', linestyle='-.', label='Raga-II', **kx)
grid[1].plot(t_raga * S2YR / 1.0e6, R_SB / PC2CM, color='purple', linestyle='--', dashes=(4,4), label='STARBENCH', **kx2)

grid[0].plot(t * S2YR / 1.0e6, error_SB, color='purple', label='SB', **kx2)

### Legend
handles, labels = grid[1].get_legend_handles_labels()
legend = grid[0].legend(handles, labels, loc=(0.05, 0.45), fontsize=10)
legend.get_frame().set_linewidth(plotter.linewidth)

###	Save figure.
with warnings.catch_warnings():
	warnings.simplefilter("ignore")
	plotter.save_plot(outputfile)

print sys.argv[0] + ': plotted ' + outputfile

