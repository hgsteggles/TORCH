import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy.io import fits
import math

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData(useSet4=True)
masses = mp1_data.masses

load_src("koo", "../torchpack/koo.py")
import koo

YR2S = 3.154e7
S2YR = 1.0 / YR2S
PC2CM = 3.086e18
CM2PC = 1.0 / PC2CM

params = np.genfromtxt("config/zams-stripped.txt", skip_header=1)

massdat = params[:,0] # [Msun]
logQdat = params[:,6] # log[Q / phot.s-1]

logQ_interp = interpolate.interp1d(massdat, logQdat)

nhs = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]

def calcStagnationTime(logQ, nh):
	ci = koo.soundSpeed(8000.0, 0.5)
	co = koo.soundSpeed(10.0, 1.0)
	cico = ci / co

	RS = koo.calcStromgrenRadius(logQ, nh)
	ts = RS / ci

	return (4.0 / 7.0) * ((3.0 / 4.0)**0.5) * (((4.0 / 3.0)**(7.0 / 6.0)) * (cico**(7.0 / 3.0)) - 1.0) * ts

for i in range(45):
	imass = i % 9
	inh = int(i / 9)

	nh = nhs[inh]
	logQ = logQ_interp(masses[imass])

	t_stag = calcStagnationTime(logQ, nh) * S2YR / 1000.0

	print str(i) + " " + str(t_stag)



