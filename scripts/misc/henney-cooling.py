import math
import numpy as np
import scipy.interpolate
import scipy
import matplotlib.pyplot as plt

KB = 1.3806485e-16

def extrap1d(interpolator, x):
	xs = interpolator.x
	ys = interpolator.y

	if x < xs[0]:
		return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
	elif x > xs[-1]:
		return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
	else:
		return interpolator(x)

def coolCollExcIonMetals(nhii, T):
	T1 = 33610.0
	T2 = 2180.0
	z0 = 5.0e-4
	return 2.905e-19 * z0 * (nhii**2) * math.exp(-((T1 / T) + (T2 / T)**2))

def coolCollExcNeuMetals(nhii, nhi, T):
	T3 = 28390.0
	T4 = 1780.0
	z0 = 5.0e-4
	return 4.477e-20 * z0 * nhii * nhi * math.exp(-((T3 / T) + (T4 / T)**2))

def getHummerSpline():
	data = np.genfromtxt("refdata/hummer.txt",skip_header=1)
	logT = data[:,0]
	T = np.power(10.0, logT)

	bB = data[:,4] / np.sqrt(T)
	return scipy.interpolate.interp1d(T, bB, kind='linear')

def getHummerSpline2():
	data = np.genfromtxt("refdata/hummer.txt",skip_header=1)
	logT = data[:,0]
	bBsqrtT = data[:,6]

	return scipy.interpolate.interp1d(logT, bBsqrtT, kind='linear')

def getRagaSpline():
	T = [3162.2776602, 3981.0717055, 5011.8723363, 6309.5734448, 7943.2823472,
		10000.0000000, 12589.2541179, 15848.9319246, 19952.6231497, 25118.8643151,
		31622.7766017, 39810.7170553, 50118.7233627, 63095.7344480, 79432.8234724,
		100000.0000000, 125892.5411794, 158489.3192461, 199526.2314969, 251188.6431510,
		316227.7660168, 398107.1705535, 501187.2336273, 630957.3444802, 794328.2347243,
		1000000.0000000]
	R = [1.150800e-34, 2.312065e-31, 9.571941e-29, 1.132400e-26,	4.954502e-25,
		9.794900e-24, 1.035142e-22, 6.652732e-22, 2.870781e-21, 9.036495e-21, 2.218196e-20,
		4.456562e-20, 7.655966e-20, 1.158777e-19, 1.588547e-19, 2.013724e-19, 2.393316e-19,
		2.710192e-19, 2.944422e-19, 3.104560e-19, 3.191538e-19, 3.213661e-19, 3.191538e-19,
		3.126079e-19, 3.033891e-19, 2.917427e-19]
	return scipy.interpolate.interp1d(np.log10(T), np.log10(R), kind='linear')

ragaSpline = getRagaSpline()
hummerSpline = getHummerSpline2()

def coolFreeFree(nhii, T):
	#2.59e-13 * pow(T / 10000.0, -0.7)
	rate = extrap1d(hummerSpline, np.log10(T)) / np.sqrt(T)
	return (nhii**2) * KB * T * rate

def coolCollExcNeuHydrogen(nhii, nhi, T):
	damp = 5.0e5
	rate = extrap1d(ragaSpline, math.log10(T))
	return nhi * nhii * math.exp((2.302585093 * rate) - (T / damp)**2)

def coolCIE(nhii, T):
	minT = 5.0e4
	z0 = 5.0e-4
	if T > minT:
		rate = 3.485e-15 * nhii**2 * z0 * math.exp(-0.63 * math.log(T)) * (1.0 - math.exp(-math.pow(1.0e-5*T, 1.63)))
		smoothing = min(1.0, (T - 5.0e4) / (2.0e4))
		return rate * smoothing
	else:
		return 0

def coolPDR(nh, T):
	T0 = 70.0 + 220.0 * math.pow(nh / 1.0e6, 0.2)
	return 3.981e-27 * (nh**1.6) * (T**0.5) * math.exp(-T0 / T)

def coolTotal(hii, nh, T):
	nhii = hii * nh
	nhi = (1.0 - hii) * nh

	mp = coolCollExcIonMetals(nhii, T)
	mn = coolCollExcNeuMetals(nhii, nhi, T)
	hp = coolFreeFree(nhii, T)
	hn = coolCollExcNeuHydrogen(nhii, nhi, T)
	cie = coolCIE(nhii, T)
	pdr = coolPDR(nh, T)

	return mp + mn + hn + hp + cie + pdr

def plot():
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
	outputfile_qlyc = "henney.png"
	plot_type = "if" # ["if", "hii", "ssw"]

	torch.set_font_sizes(fontsize)

	### Data
	lmp = []
	lmn = []
	lhp = []
	lhn = []
	lcie = []
	lpdr = []
	ltot = []

	if plot_type == "ssw":
		nh = 1
	elif plot_type == "hii":
		nh = 1000
	else:
		nh = 10000

	hii = 1.0
	if plot_type == "if":
		hii = 0.5
	N = 1001
	minlogT = 4.0
	maxlogT = 10.0
	if plot_type == "if" or plot_type == "hii":
		minlogT = 3.8
		maxlogT = 4.0
	isLog = True

	nhii = hii * nh
	nhi = (1.0 - hii) * nh

	if not isLog:
		tem = np.linspace(10.0**minlogT, 10.0**maxlogT, N, endpoint=True)
	else:
		tem = np.linspace(minlogT, maxlogT, N, endpoint=True)

	for i in range(len(tem)):
		itemp = tem[i]
		if isLog:
			itemp = math.pow(10.0, itemp)
		lmp.append(coolCollExcIonMetals(nhii, itemp))
		lmn.append(coolCollExcNeuMetals(nhii, nhi, itemp))
		lhp.append(coolFreeFree(nhii, itemp))
		lhn.append(coolCollExcNeuHydrogen(nhii, nhi, itemp))
		lcie.append(coolCIE(nhii, itemp))
		lpdr.append(coolPDR(nhi, itemp))
		ltot.append(lmp[i] + lmn[i] + lhp[i] + lhn[i] + lcie[i] + lpdr[i])

	### Plotting.
	plotter = torch.Plotter(1, 1, plot_size, figformat, DPI)

	###	Axes.
	grid = plotter.axes1D((1,1), aspect_ratio=0.75)
	grid[0].set_xlabel(plotter.format_label(torch.VarType('T', units='K', isLog10=True)))
	grid[0].set_ylabel(plotter.format_label(torch.VarType('L', units='erg\\:s^{-1}\\:cm^{-3}', isLog10=False)))

	### Plot.
	if plot_type == "if":
		grid[0].plot(tem, lmn, label=r"$L_\mathrm{M^0}$", color='m')
		grid[0].plot(tem, lhn, label=r"$L_\mathrm{H^0}$", color='y')
		grid[0].plot(tem, lmp, label=r"$L_\mathrm{M^+}$", color='g')
		grid[0].plot(tem, lhp, label=r"$L_\mathrm{H^+}$", color='r')
		grid[0].plot(tem, ltot, label=r"$L_\mathrm{tot}$", color='k')
	elif plot_type == "hii":
		grid[0].plot(tem, lmp, label=r"$L_\mathrm{M^+}$", color='g')
		grid[0].plot(tem, lhp, label=r"$L_\mathrm{H^+}$", color='r')
		grid[0].plot(tem, ltot, label=r"$L_\mathrm{tot}$", color='k')
	elif plot_type == "ssw":
		grid[0].plot(tem, lmp, label=r"$L_\mathrm{M^+}$", color='g')
		grid[0].plot(tem, lhp, label=r"$L_\mathrm{H^+}$", color='r')
		grid[0].plot(tem, lcie, label=r"$L_\mathrm{CIE}$", color='b')
		grid[0].plot(tem, ltot, label=r"$L_\mathrm{tot}$", color='k')

	grid[0].set_xlim([minlogT, maxlogT])
	if plot_type == "ssw":
		grid[0].legend(loc='upper right')
	else:
		grid[0].legend(loc='upper left')

	###	Save figure.
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		plotter.save_plot(outputfile_qlyc)

	print sys.argv[0] + ': plotted in ' + outputfile_qlyc

plot()