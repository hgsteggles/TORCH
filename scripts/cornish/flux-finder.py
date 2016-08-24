import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
from astropy.io import fits
import math

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Finds closest matching parameters for flux.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("koo", "../torchpack/koo.py")
import koo

YR2S = 3.154e7
S2YR = 1.0 / YR2S
S2KYR = S2YR / 1000.0
PC2CM = 3.086e18
CM2PC = 1.0 / PC2CM

params = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

massdat = params[:,0] # [Msun]
logQdat = params[:,6] # log[Q / phot.s-1]

logQ_interp = interpolate.interp1d(massdat, logQdat)

def calcPlanck2(tem, freq):
	l = 3.0e10 / freq
	c1 = 3.7417749e-5
	c2 = 1.4387687
	val = c2 / (l * tem)
	return (c1 / (l**5)) / (math.exp(val) - 1) * 1.0e-11

def calcTau(T, freq, ne, logq):
	return 0.53 * (T / 1.0e4)**(-1.35) * (freq / 1.0e9)**(-2.1) * (ne / 1.0e3)**0.66 * 10.0**((logq - 49.7) / 3.0)

nhs = np.arange(0.8e4, 20.0e4, 0.1e4)
masses = np.arange(6.0, 70.0, 1.0)
ages = np.arange(2.0, 202.0, 2.0) * 1000.0 * YR2S
dsun = 14.0 * 1000.0 * PC2CM

T = 8000.0
f = 5.0e9
B = calcPlanck2(T, f)

targetF = 100.0
diffmin = 1.0e30
jnh = 0
jmass = 0
jage = 0

for imass in range(len(masses)):
	for iage in range(len(ages)):
		for inh in range(len(nhs)):
			age = ages[iage]
			mass = masses[imass]
			nh = nhs[inh]
			logQ = logQ_interp(mass)

			stro = koo.calcStromgrenRadius(logQ, nh)
			raga = koo.calcSpitzerRadius(logQ, nh, age)
			stag = koo.calcStagnationRadius(logQ, nh)
			radius = np.minimum(raga, stag)

			angsize_rad = 2.0 * (radius / dsun)
			angsize_as = angsize_rad * 60.0 * 60.0 * 180.0 / math.pi
			ne = nh * (stro / radius)**(1.5)
			tau = calcTau(T, f, ne, logQ)
			S = B * (1.0 - math.exp(-tau))

			F = S * (0.5 * angsize_rad)**2 * 1.0e26 * 1000

			diff = abs(F - targetF)

			if diff < diffmin:
				jmass = imass
				jage = iage
				jnh = inh
				diffmin = diff


print "nh  = " + str(nhs[jnh])
print "M   = " + str(masses[jmass])
print "Age = " + str(ages[jage] * S2KYR)