import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy.io import fits
import math

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
parser.add_argument('-ben', action='store_true', help='Use simple ben stromgren survey.')
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
nstars = len(star_data[:,0])

params = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

massdat = params[:,0] # [Msun]
logQdat = params[:,6] # log[Q / phot.s-1]

logQ_interp = interpolate.interp1d(massdat, logQdat)

benstr = ""
if args.ben:
	benstr = "-ben"

### Data
simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + benstr + ".txt", skip_header=1)

survey_names = np.genfromtxt("data/cornish/cornish-uchiis.txt", skip_header=1, delimiter=',', dtype=str)[:,0]
distance_names = np.genfromtxt("data/cornish/cornish-distances.txt", skip_header=1, dtype=str)[:,0]

phy_sizes = simulated_survey[:,7] * simulated_survey[:,4] * 1000.0 * math.pi / (60.0 * 60.0 * 180.0)

nhs = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]
nes = []

def calc_ne(mass, size):
	nh = nhs[args.iden]
	logQ = logQ_interp(mass)
	stro = koo.calcStromgrenRadius(logQ, nh)
	radius = 0.5 * size * PC2CM
	ne = nh * (stro / radius)**(1.5)

	return ne

for i in range(len(simulated_survey[:,0])):
	nes.append(calc_ne(simulated_survey[i, 1], phy_sizes[i]))

print "Simulated mean ne = " + str(np.mean(nes)) + " cm^-3"