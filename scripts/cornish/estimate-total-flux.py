import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy.io import fits
import math

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData(useSet4=True)
masses = mp1_data.masses

load_src("cdat", "../torchpack/cornishdata.py")
import cdat
cornish_data = cdat.CornishData(args.iden)

load_src("koo", "../torchpack/koo.py")
import koo

YR2S = 3.154e7
PC2CM = 3.086e18
CM2PC = 1.0 / PC2CM

star_data = cornish_data.star_data

params = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

massdat = params[:,0] # [Msun]
logQdat = params[:,6] # log[Q / phot.s-1]

logQ_interp = interpolate.interp1d(massdat, logQdat)

ofile = open(cornish_data.dirname + "/estimate-total-fluxes.txt", 'w')
ofile.write("star_id total_flux[mJy]\n")

def calcPlanck(tem, freq):
	h = 6.626e-34
	c = 3.0e8
	kb = 1.38e-23
	return (2.0 * h * freq * freq * freq / (c * c)) / (math.exp(h * freq / (kb * tem)) - 1)

def calcTau(T, freq, ne, logq):
	return 0.53 * (T / 1.0e4)**(-1.35) * (freq / 1.0e9)**(-2.1) * (ne / 1.0e3)**0.66 * 10.0**((logq - 49.7) / 3.0)

n = 0

for i in range(1, len(star_data[:,0]) + 1):
	ipad = "%03d" % (i,)
	dirname = cornish_data.dirname + "/star_" + ipad

	age = star_data[i - 1, 3] * YR2S * 1000.0

	mass = star_data[i - 1, 0]
	nh = 3.e4
	logQ = logQ_interp(mass)
	dsun = star_data[i - 1, 4] * 1000.0 * PC2CM

	stro = koo.calcStromgrenRadius(logQ, nh)
	raga = koo.calcSpitzerRadius(logQ, nh, age)
	stag = koo.calcStagnationRadius(logQ, nh)
	radius = np.minimum(raga, stag)

	angsize_rad = (radius / dsun)
	angsize_as = angsize_rad * 60.0 * 60.0 * 180.0 / math.pi
	ne = nh * (stro / radius)**(1.5)
	T = 8000.0
	f = 5.0e9
	tau = calcTau(T, f, ne, logQ)
	B = calcPlanck(T, f)
	S = B * (1.0 - math.exp(-tau))

	F = S * math.pi * angsize_rad * angsize_rad * 1.0e26 * 1000

	if F > 1000 and angsize_as < 24:
		n += 1
		print str(i) + " " + str(F) + " " + str(angsize_as)

print n

ofile.close()