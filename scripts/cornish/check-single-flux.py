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
PC2CM = 3.086e18
CM2PC = 1.0 / PC2CM

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

age = 80 * 1000.0 * YR2S
mass = 10.67
nh = 6.4e4
logQ = logQ_interp(mass)
dsun = 4.8251663 * 1000.0 * PC2CM

stro = koo.calcStromgrenRadius(logQ, nh)
raga = koo.calcSpitzerRadius(logQ, nh, age)
stag = koo.calcStagnationRadius(logQ, nh)
radius = np.minimum(raga, stag)

angsize_rad = 2.0 * (radius / dsun)
angsize_as = angsize_rad * 60.0 * 60.0 * 180.0 / math.pi
ne = nh * (stro / radius)**(1.5)

T = 8000.0
f = 5.0e9
tau = calcTau(T, f, ne, logQ)

B1 = calcPlanck(T, f) * 1.0e26
B2 = calcPlanck2(T, f) * 1.0e26

print str(B1) + " " + str(B2)

S = B2 * (1.0 - math.exp(-tau))

F = S * (0.5 * angsize_rad)**2 * 1000

print str(F) + " " + str(angsize_as)