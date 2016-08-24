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

def calcPlanck(tem, freq):
	h = 6.626e-34
	c = 3.0e8
	kb = 1.38e-23
	return (2.0 * h * freq * freq * freq / (c * c)) / (math.exp(h * freq / (kb * tem)) - 1)

def calcTau(T, freq, ne, logq):
	return 0.53 * (T / 1.0e4)**(-1.35) * (freq / 1.0e9)**(-2.1) * (ne / 1.0e3)**0.66 * 10.0**((logq - 49.7) / 3.0)

sizes = np.genfromtxt(cornish_data.dirname + "/sizes.txt", skip_header=1)

nhs = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]

for i in range(1, len(star_data[:,0]) + 1):
	ipad = cornish_data.star_id_str(i)
	dirname = cornish_data.dirname + "/star_" + ipad + "/"

	age = star_data[i - 1, 3] * 1000.0 * YR2S
	mass = star_data[i - 1, 0]
	nh = nhs[args.iden]
	logQ = logQ_interp(mass)
	dsun = star_data[i - 1, 4] * 1000.0 * PC2CM

	stro = koo.calcStromgrenRadius(logQ, nh)
	raga = koo.calcSpitzerRadius(logQ, nh, age)
	stag = koo.calcStagnationRadius(logQ, nh)
	radius = np.minimum(raga, stag)

	angsize_rad = (radius / dsun)
	angsize_as = 2.0 * angsize_rad * 60.0 * 60.0 * 180.0 / math.pi

	if 0.75 * angsize_as > sizes[i - 1][1]:
		print i