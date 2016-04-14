import numpy as np
import math


CM2PC = 3.24078e-19
PC2CM = 1.0/CM2PC
R_GAS = 8.3144598e7
mh = 1.6738232e-24
S2YR = 3.17098e-8
YR2S = 3.154e7

def calcStromgren(logQ, nh):
	T = 8000
	alphaB = 2.59e-13*math.pow(T/10000.0, -0.7)
	return pow(3.0 * pow(10.0, logQ) / (4.0 * math.pi * nh * nh * alphaB), 1.0/3.0)

def calcPressureEqual2(mdot, vinf, rho_0):
	T = 8000.0
	c5 = math.pow(math.sqrt(R_GAS * T / 0.5), 5.0)
	return 0.142 * math.sqrt(0.5 * mdot * vinf * vinf / (rho_0 * c5))

def finalRadius(mdot, vinf, rho_0):
	tp25 = math.pow(calcPressureEqual2(mdot, vinf, rho_0), 2.0/5.0)
	return 0.92 * math.pow(0.5 * mdot * vinf * vinf / (rho_0 * math.pow(vinf, 5.0/3.0)), 0.3) * tp25

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

for i in range(len(den)):
	for j in range(len(mass)):
		id = "{:2d}".format(9*i + j + 1)
		m = "{:3d}".format(int(mass[j]))
		t = "{:1.2f}".format(temp[j])
		r = "{:4.2f}".format(rstar[j])
		l = "{:4.2f}".format(logL[j])
		q = "{:4.2f}".format(logQ[j])
		md = "{:8.2e}".format(mdot[j] / 6.341958397e25) # [Msun.yr-1]
		v = "{:8.2e}".format(vinf[j] / 100000.0) # [km.s-1]
		d = "{:8.2e}".format(den[i])
		ysc = "{:3.2f}".format(calcStromgren(logQ[j], den[i]) / cloud_dist)
		tp = "{:5.2f}".format(calcPressureEqual2(mdot[j], vinf[j], mh * den[i]) * S2YR / 1000.0)
		rfinal = "{:4.3f}".format(finalRadius(mdot[j], vinf[j], mh * den[i]) * CM2PC)

		print id + " & " + m + " & " + t + " & " + r + " & " + l + " & " + q + " & " + md \
			  + " & " + v + " & " + d + " & " + ysc + " & " + tp + " & " + rfinal + " \\\\"

