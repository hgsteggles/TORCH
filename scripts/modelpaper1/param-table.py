import numpy as np
import math
import sys

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("koo", "../torchpack/koo.py")
load_src("henney", "../misc/henney-cooling.py")

import koo
import henney

import lupa
lua = lupa.LuaRuntime(unpack_returned_tuples=True)

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData()

CM2PC = 3.24078e-19
PC2CM = 1.0/CM2PC
R_GAS = 8.3144598e7
mh = 1.6738232e-24
S2YR = 3.17098e-8
YR2S = 3.154e7
KB = 1.3806485e-16

def critL(pre0, nh0, r_inj, vinf):
	g = 5.0 / 3.0
	#T = 147910.0
	mu = 0.5
	#T = 1.0e8
	#c = math.sqrt(R_GAS * T / mu)
	c = math.sqrt(g * pre0 / (mh * nh0))
	T = pre0 * mu / (R_GAS * mh * nh0)

	alpha = 0.28
	#COOL = 3.185e-15 * math.pow(T, -0.63) * (1.0 - math.exp(-math.pow(T / 1.0e5, 1.63)))
	COOL = henney.coolTotal(1.0, nh0, T)
	coeff1 = 6.0 * (g - 1.0) * math.pi * alpha**2 * mu**2 * mh**2 * r_inj * vinf**4
	coeff2 = (g + 1.0) * COOL
	return (coeff1 / coeff2) * ((vinf**2 / 2.0) - (c**2 / (g - 1.0)))

def critL2(r_inj, vinf):
	return 4.4e40 * (r_inj * CM2PC) * (vinf / 1.0e8)

def coolTime(pre0, nh0):
	mu = 0.5
	T = 1.0e8
	T = pre0 * mu / (R_GAS * mh * nh0)
	COOL = henney.coolTotal(1.0, nh0, T)
	return 1.5 * nh0 * KB * T / (COOL)

def exitTimeCoeff(r_inj, vinf):
	g = 5.0 / 3.0
	A = math.sqrt((g - 1) / (g + 1)) * math.pow((g + 1) / (6.0 * g + 2), (3.0 * g + 1) / (5.0 * g + 1))
	return (r_inj / vinf) / A

def exitTime(r_inj, vinf):
	g = 5.0 / 3.0
	cinj = math.sqrt((g - 1) / (g + 1)) * vinf
	return r_inj / cinj

def initialRadius(pre0, nh0, rinj, vinf):
	return - coolTime(pre0, nh0) / exitTimeCoeff(r_inj, vinf)


fileprefix2 = "data/model_paper1/set2/"
fileprefix3 = "data/model_paper1/set3/"

snapshots = [10, 20, 30, 40, 50]

nxs = []
nys = []
sizexs = []

for i in range(45):
	col = int(i / 9)
	row = i % 9
	itime = 50 if i >= 19 and i <= 27 else 20

	filename = mp1_data.getParamFilename(row, col, itime)
	filestring = open(filename, 'r').read()

	table = lua.eval("{" + filestring + "}")
	nxs.append(table["Parameters"]["Grid"]["no_cells_x"])
	nys.append(table["Parameters"]["Grid"]["no_cells_y"])
	sizexs.append(table["Parameters"]["Grid"]["side_length"])

###	Data set up.
data = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

mass = data[:,0] # [Msun]
temp = data[:,1]
logQ = data[:,6] # log[Q / phot.s-1]
mdot = data[:,8] # [g.s-1]
vinf = data[:,9] # [cm.s-1]
logL = data[:,2] # log[L / Lsun]
rstar = data[:,3] # [Rsun]

den = np.array([0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4])

cloud_dist = 0.35 * PC2CM

for j in range(len(mass)):
	m = "{:3d}".format(int(mass[j]))
	t = "{:5.0f}".format(temp[j])
	r = "{:4.2f}".format(rstar[j])
	l = "{:4.2f}".format(logL[j])
	q = "{:4.2f}".format(logQ[j])
	ed = "{:4.2e}".format(0.5 * mdot[j] * vinf[j] * vinf[j]) # erg.s-1
	md = "{:8.2e}".format(mdot[j] / 6.341958397e25) # [Msun.yr-1]
	v = "{:8.2e}".format(vinf[j] / 100000.0) # [km.s-1]

	print m + " & " + t + " & " + r + " & " + l + " & " + q + " & " + ed + " & " + md \
			  + " & " + v + " \\\\"

for i in range(len(den)):
	for j in range(len(mass)):
		index = 9*i + j
		m = "{:3d}".format(int(mass[j]))
		d = "{:8.2e}".format(den[i])
		ysc = "{:4.3f}".format(koo.calcStromgrenRadius(logQ[j], den[i]) / cloud_dist)
		nx = str(nxs[index])
		ny = str(nys[index])
		nxny = "{($" + nx + "$,$" + ny + "$)}"
		sizex = "{:3.2f}".format(sizexs[index] / 3.09e18)
		r_inj = 10.0 * sizexs[index] / nxs[index]
		lcrit = "{:4.3e}".format(critL2(r_inj, vinf[j]))
		tstart = "{:4.3f}".format(koo.calcStartTime(mdot[j], vinf[j], mh*den[i], 8000.0, 0.5, r_inj) * S2YR / 1000.0)
		tcoolval = coolTime(mp1_data.central_pressures[index], mp1_data.central_densities[index]) * S2YR / 1000.0
		tcool = "{:4.3e}".format(tcoolval)
		texitval = exitTime(r_inj, vinf[j]) * S2YR / 1000.0
		texit = "{:4.3f}".format(texitval)
		R_inj = "{:4.3f}".format(r_inj * CM2PC)
		ri = "{:8d}".format((int)(-tcoolval / texitval))

		print m + " & " + d + " & " + nx + " & " + ny + " & " + sizex + " & " + ysc + " & " + R_inj + " & " + lcrit + " & " + tstart + " & " + tcool + " & " + ri + " \\\\"