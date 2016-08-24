import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.io import fits
import math
from scipy import interpolate

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Calculates accretion luminosity of stars.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("cdat", "../torchpack/cornishdata.py")
import cdat
cornish_data = cdat.CornishData(args.iden)
star_data = cornish_data.star_data
nstars = len(star_data[:,0])

params = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

massdat = params[:,0] # [Msun]
Rdat = params[:,3]  # [Rsun]

Rdat_interp = interpolate.interp1d(massdat, Rdat)

mcur = star_data[:,0]
mfin = star_data[:,1]
t_ms = star_data[:,2]
lbol = star_data[:,10]

for i in range(nstars):
	if mcur[i] > 20 and mcur[i] < mfin[i]:
		G = 1.0e-7
		SIG_cl = 1.0
		MSUN = 2e33
		YEAR = 3.154e7
		RSUN = 7.0e10
		LSUN = 3.828e33

		mt_acc = 4.6e-4 * math.pow(mfin[i] / 30.0, 0.75) * math.pow(SIG_cl, 0.75) * math.pow(mcur[i] / mfin[i], 0.5) * MSUN / YEAR

		L_acc = G * mcur[i] * MSUN * mt_acc / (Rdat_interp(mcur[i]) * RSUN)

		print math.log10(L_acc / LSUN)