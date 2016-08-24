import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import multiprocessing
import signal
from scipy import interpolate

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData(useSet4=True)
masses = mp1_data.masses
nhs = mp1_data.densities
injection_radii = mp1_data.injection_radii

params = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

massdat = params[:,0] # [Msun]
logQdat = params[:,6] # log[Q / phot.s-1]

logQ_interp = interpolate.interp1d(massdat, logQdat)

KB = 1.38e-16
RGAS = 8.3144598e7

YR2S = 3.154e7
S2YR = 1.0 / YR2S
PC2CM = 3.086e18
CM2PC = 1.0 / PC2CM

def soundSpeed(T, mu):
	return math.sqrt(RGAS * T / mu)

def calcAlphaB(T):
	return 2.59e-13*math.pow(T / 10000.0, -0.7)

def calcStromgrenRadius(logQ, nh):
	alphaB = calcAlphaB(8000.0)
	return pow(3.0 * pow(10.0, logQ) / (4.0 * math.pi * nh * nh * alphaB), 1.0/3.0)

def calcStagnationRadius(logQ, nh):
	T = 8000.0
	ci = soundSpeed(8000.0, 0.5)
	co = soundSpeed(300.0, 1.0)
	cico = ci / co
	return math.pow(math.sqrt(4.0 / 3.0) * cico, 4.0 / 3.0) * calcStromgrenRadius(logQ, nh)

print "index inj_rad_time[kyr] inj_rad[pc] stag_rad[pc]"

for i in range(45):
	imass = i % 9
	inh = int(i / 9)

	inj_rad = injection_radii[i] * PC2CM
	mass = masses[imass]
	nh = nhs[inh]
	logQ = logQ_interp(mass)

	RS = calcStromgrenRadius(logQ, nh)
	cs = soundSpeed(8000.0, 0.5)
	ts = RS / cs
	t = (4.0 / 7.0) * np.sqrt(3.0 / 4.0) * (pow(inj_rad / RS, 7.0 / 4.0) - 1.0) * ts

	print str(i + 1) + " " + str(t * S2YR / 1000.0) + " " + str(injection_radii[i]) + " " + str(calcStagnationRadius(logQ, nh) * CM2PC)

