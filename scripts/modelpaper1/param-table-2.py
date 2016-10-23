import numpy as np
import math

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("koo", "../torchpack/koo.py")
import koo

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData()

CM2PC = 3.24078e-19
PC2CM = 1.0/CM2PC
R_GAS = 8.3144598e7
mh = 1.6738232e-24
S2YR = 3.17098e-8
YR2S = 3.154e7

mvd_iden = []
mvd_iden.append([130.0,   150.0,  250.0,  400.0,  650.0,  900.0, 1200.0, 1800.0,  2200.0])
mvd_iden.append([280.0,   250.0,  380.0,  560.0, 1000.0, 1300.0, 1800.0, 2000.0,  2200.0])
mvd_iden.append([610.0,   480.0,  650.0,  900.0, 1300.0, 1800.0, 2200.0, 2400.0,  2700.0])
mvd_iden.append([800.0,   900.0, 1000.0, 1400.0, 2100.0, 2500.0, 3000.0, 5000.0,  5000.0])
mvd_iden.append([2000.0, 2200.0, 2200.0, 2400.0, 3000.0, 4500.0, 5500.0, 7000.0, 12000.0])

mvt_iden = []
mvt_iden.append([520.0, 620.0, 1000.0, 1700.0, 2600.0, 4500.0, 6500.0, 9000.0, 12000.0])
mvt_iden.append([620.0, 520.0,  750.0, 1000.0, 1600.0, 2400.0, 2500.0, 5000.0,  8000.0])
mvt_iden.append([580.0, 520.0,  600.0,  800.0, 1400.0, 1800.0, 1900.0, 2500.0,  3000.0])
mvt_iden.append([520.0, 550.0,  550.0,  700.0, 1100.0, 1300.0, 1500.0, 1700.0,  1800.0])
mvt_iden.append([520.0, 540.0,  520.0,  620.0,  850.0, 1200.0, 1300.0, 1400.0,  1500.0])

###	Data set up.
data = np.genfromtxt("config/zams-stripped.txt", skip_header=1)

mass = data[:,0] # [Msun]
temp = data[:,1]
logQ = data[:,6] # log[Q / phot.s-1]
mdot = data[:,8] # [g.s-1]
vinf = data[:,9] # [cm.s-1]
logL = data[:,2] # log[L / Lsun]
rstar = data[:,3] # [Rsun]

den = mp1_data.densities

cloud_dist = 0.35 * PC2CM

for i in range(len(den)):
	for j in range(len(mass)):
		index = 9*i + j + 1
		id = "{:2d}".format(index)
		lower = "{:3.2f}".format(koo.calcLowerPRB(vinf[j]))

		vin_vcr = "{:6.3f}".format(vinf[j] / koo.calcCritVel(mdot[j], vinf[j], mh * mvd_iden[i][0]))
		tp_rb = "{:5.2f}".format(koo.calcConfineTimeRB(mdot[j], vinf[j], mh * mvd_iden[i][0], 8000.0, 0.5) * S2YR / 1000.0)
		rfinal_rb = "{:4.3f}".format(koo.calcFinalRadiusRB(mdot[j], vinf[j], mh * mvd_iden[i][0], 8000.0, 0.5) * CM2PC)
		tp_prb = "{:5.2f}".format(koo.calcConfineTimePRB(mdot[j], vinf[j], mh * mvd_iden[i][0], 8000.0, 0.5) * S2YR / 1000.0)
		rfinal_prb = "{:4.3f}".format(koo.calcFinalRadiusPRB(mdot[j], vinf[j], mh * mvd_iden[i][0], 8000.0, 0.5) * CM2PC)
		tp_ab = "{:5.2f}".format(koo.calcConfineTimeAB(mdot[j], vinf[j], mh * mvd_iden[i][0], 8000.0, 0.5) * S2YR / 1000.0)
		rfinal_fab = "{:4.3f}".format(koo.calcFinalRadiusFAB(mdot[j], vinf[j], mh * mvd_iden[i][0], 8000.0, 0.5) * CM2PC)
		rfinal_pab = "{:4.3f}".format(koo.calcFinalRadiusPAB(mdot[j], vinf[j], mh * mvd_iden[i][0], 8000.0, 0.5) * CM2PC)

		vin_vcr2 = "{:6.3f}".format(vinf[j] / koo.calcCritVel(mdot[j], vinf[j], mh * den[i]))
		tp_rb2 = "{:5.2f}".format(koo.calcConfineTimeRB(mdot[j], vinf[j], mh * den[i], 300.0, 1.0) * S2YR / 1000.0)
		rfinal_rb2 = "{:4.3f}".format(koo.calcFinalRadiusRB(mdot[j], vinf[j], mh * den[i], 300.0, 1.0) * CM2PC)
		tp_prb2 = "{:5.2f}".format(koo.calcConfineTimePRB(mdot[j], vinf[j], mh * den[i], 300.0, 1.0) * S2YR / 1000.0)
		rfinal_prb2 = "{:4.3f}".format(koo.calcFinalRadiusPRB(mdot[j], vinf[j], mh * den[i], 300.0, 1.0) * CM2PC)
		tp_ab2 = "{:5.2f}".format(koo.calcConfineTimeAB(mdot[j], vinf[j], mh * den[i], 300.0, 1.0) * S2YR / 1000.0)
		rfinal_fab2 = "{:4.3f}".format(koo.calcFinalRadiusFAB(mdot[j], vinf[j], mh * den[i], 300.0, 1.0) * CM2PC)
		rfinal_pab2 = "{:4.3f}".format(koo.calcFinalRadiusPAB(mdot[j], vinf[j], mh * den[i], 300.0, 1.0) * CM2PC)

		t_start = str(koo.calcStartTime(mdot[j], vinf[j], mh * den[i], 8000.0, 0.5, mp1_data.injection_radii[index-1] * PC2CM) * S2YR / 1000.0)

		print id + " & " + lower + " & " +\
			  vin_vcr + " & " + tp_rb + " & " + rfinal_rb + " & " + tp_prb + " & " +\
			  rfinal_prb + " & " + tp_ab + " & " + rfinal_fab + " & " + rfinal_pab + " & " +\
			  vin_vcr2 + " & " + tp_rb2 + " & " + rfinal_rb2 + " & " + tp_prb2 + " & " +\
			  rfinal_prb2 + " & " + tp_ab2 + " & " + rfinal_fab2 + " & " + rfinal_pab2 + " & " +\
			  t_start + "\\\\"

