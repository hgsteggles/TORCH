import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import lupa
lua = lupa.LuaRuntime(unpack_returned_tuples=True)

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
nstars = len(star_data[:,0])

params = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

massdat = params[:,0] # [Msun]
logQdat = params[:,6] # log[Q / phot.s-1]

logQ_interp = interpolate.interp1d(massdat, logQdat)

def calcPlanck(tem, freq):
	h = 6.626e-34
	c = 3.0e8
	kb = 1.38e-23
	return (2.0 * h * freq * freq * freq / (c * c)) / (math.exp(h * freq / (kb * tem)) - 1)

def calcPlanck2(tem, freq):
	l = 3.0e10 / freq
	c1 = 3.7417749e-5
	c2 = 1.4387687
	val = c2 / (l * tem)
	return (c1 / (l**5)) / (math.exp(val) - 1) * 1.0e-11

def calcTau(T, freq, ne, logq):
	return 0.53 * (T / 1.0e4)**(-1.35) * (freq / 1.0e9)**(-2.1) * (ne / 1.0e3)**0.66 * 10.0**((logq - 49.7) / 3.0)

nhs = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]

mcur = star_data[:,0]
t_ms = star_data[:,2]
age = star_data[:,3]
d_sun = star_data[:,4]
phy_size = np.zeros(nstars, dtype=float)
ang_size = np.zeros(nstars, dtype=float)
casa_ang_size = np.zeros(nstars, dtype=float)
tot_flux = np.zeros(nstars, dtype=float)
casa_tot_flux = np.zeros(nstars, dtype=float)
casa_peak_flux = np.zeros(nstars, dtype=float)
g_long = star_data[:,15]
g_lat = star_data[:,16]

for i in range(1, nstars + 1):
	t = age[i - 1] * 1000.0 * YR2S
	m = mcur[i - 1]
	nh = nhs[args.iden]
	logQ = float(logQ_interp(m))
	dsun = d_sun[i - 1] * 1000.0 * PC2CM

	stro = koo.calcStromgrenRadius(logQ, nh)
	raga = koo.calcSpitzerRadius(logQ, nh, t)
	stag = koo.calcStagnationRadius(logQ, nh)
	radius = np.minimum(raga, stag)

	angsize_rad = (2.0 * radius / dsun)
	angsize_as = angsize_rad * 60.0 * 60.0 * 180.0 / math.pi
	ne = nh * (stro / radius)**(1.5)
	T = 8000.0
	f = 5.0e9
	tau = calcTau(T, f, ne, logQ)
	B = calcPlanck2(T, f) * 1.0e26
	S = B * (1.0 - math.exp(-tau))
	F = S * (0.5 * angsize_rad)**2 * 1000

	bmsize = 1.5

	pkflux = F / (max(angsize_as, bmsize) / bmsize)**2

	phy_size[i - 1] = angsize_rad * d_sun[i - 1]
	ang_size[i - 1] = angsize_as
	casa_ang_size[i - 1] = angsize_as
	tot_flux[i - 1] = F
	casa_tot_flux[i - 1] = F
	casa_peak_flux[i - 1] = pkflux

remaining_set = np.ones(nstars + 1, dtype=bool)

for i in range(1, nstars + 1):
	# Large size
	if casa_ang_size[i - 1] > 24:
		remaining_set[i - 1] = False

	# Weak flux
	ipad = cornish_data.star_id_str(i)
	dirname = cornish_data.dirname + "/star_" + ipad + "/"
	filestring = open(dirname + "radioconfig_0.lua", 'r').read()
	table = lua.eval("{" + filestring + "}")
	dec = table["Parameters"]["declination"]
	if (dec >= 14.2 and casa_peak_flux[i - 1] < 5 * 0.25) or (dec < 14.2 and casa_peak_flux[i - 1] < 5 * 0.35):
		remaining_set[i - 1] = False

ofile = open(cornish_data.dirname + "/final-survey-harry.txt", 'w')
ofile.write("star_id mcur[Msun] t_ms[kyr] age[kyr] d_sun[kpc] phy_size[pc] ang_size[\"] casa_ang_size[\"] tot_flux[mJy] casa_tot_flux[mJy] casa_pk_flux[mJy.beam-1] g_long[deg] g_lat[deg]\n")

for i in range(1, nstars + 1):
	if (remaining_set[i]):
		ofile.write(str(i) + " " + str(mcur[i-1]) + " " + str(t_ms[i-1]) +  " " + str(age[i-1])
					+ " " + str(d_sun[i-1])
					+ " " + str(phy_size[i-1]) + " " + str(ang_size[i-1])
					+ " " + str(casa_ang_size[i-1]) + " " + str(tot_flux[i-1])
					+ " " + str(casa_tot_flux[i-1]) + " " + str(casa_peak_flux[i-1])
					+ " " + str(g_long[i-1]) + " " + str(g_lat[i-1]) + '\n')

ofile.close()
