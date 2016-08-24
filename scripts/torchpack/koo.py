import math
import numpy as np

KB = 1.38e-16
RGAS = 8.3144598e7

def soundSpeed(T, mu):
	return math.sqrt(RGAS * T / mu)

def calcOldC1():
	g = 5.0 / 3.0
	l = 1.6e-19
	muh = 2.34e-24
	x = 2.3
	return (1.0 / l) * math.sqrt(8.0 * (muh * (g - 1))**5 / (KB * x * (g + 1)**8))

def calcC1():
	g = 5.0 / 3.0
	l = 3.485e-15 * 5.0e-4
	muh = 1.0e-24
	x = 1.0
	return (1.0 / l) * math.sqrt(8.0 * (muh * (g - 1))**5 / (KB * x * (g + 1)**8))

C1 = calcC1()

def calcL(mdot, vinf):
	return 0.5 * mdot * vinf * vinf

def calcCritVel(mdot, vinf, rho_0):
	L = calcL(mdot, vinf)
	return math.pow(rho_0 * L / (2.0 * math.pi * (C1**2)), 1.0 / 11.0)

def calcFiducialRadius(mdot, vinf, rho_0):
	L = calcL(mdot, vinf)
	return math.sqrt(L / (2.0 * math.pi * rho_0 * (vinf**3)))

def calcFiducialTime(mdot, vinf, rho_0):
	return calcFiducialRadius(mdot, vinf, rho_0) / vinf

def calcTransitionTime(mdot, vinf, rho_0):
	vcr = calcCritVel(mdot, vinf, rho_0)
	return math.sqrt(1.0 / 6.0) * math.pow(vcr / vinf, 11.0) * calcFiducialTime(mdot, vinf, rho_0)

def calcStartTime(mdot, vinf, rho_0, T, mu, inj_r):
	V = (4.0 * math.pi / 3.0) * inj_r * inj_r * inj_r
	g = 5.0 / 3.0
	p_amb = RGAS * T * rho_0 / mu
	L = calcL(mdot, vinf)
	return p_amb * (g - 1.0) * V / L

###############################
### Radiative Bubble
###############################

def calcBubbleRadiusRB(mdot, vinf, rho_0, time):
	return pow(3.0 * 0.5 * mdot * vinf * time * time / (math.pi * rho_0), 0.25)

def calcConfineTimeRB(mdot, vinf, rho_0, T, mu):
	L = calcL(mdot, vinf)
	cs = soundSpeed(T, mu)
	m = vinf / cs
	return math.sqrt(3.0 * L * (m**4) / (16.0 * math.pi * rho_0 * (vinf**5)))

def calcFinalRadiusRB(mdot, vinf, rho_0, T, mu):
	return calcBubbleRadiusRB(mdot, vinf, rho_0, calcConfineTimeRB(mdot, vinf, rho_0, T, mu))

###############################
### Partially Radiative Bubble
###############################

def calcLowerPRB(vinf):
	m = vinf / soundSpeed(8000.0, 0.5)
	return 0.42 / math.pow(m, 2.0 / 11.0)

def calcBubbleRadiusPRB(mdot, vinf, rho_0, time):
	L = calcL(mdot, vinf)
	num1 = 3.0 * (7**4) * C1 * L * (vinf**5) * (time**4)
	num2 = (2**11) * math.pi * (0.85**2) * (rho_0**2)
	return math.pow(num1 / num2, 1.0 / 7.0)

def calcConfineTimePRB(mdot, vinf, rho_0, T, mu):
	L = calcL(mdot, vinf)
	c0 = soundSpeed(T, mu)
	m = vinf / c0
	num1 = 24.0 * C1 * L * (m**7.0)
	num2 = (7**3) * math.pi * (0.85**2) * (rho_0**2.0) * (vinf**2)
	return math.pow(num1 / num2, 1.0 / 3.0)

def calcFinalRadiusPRB(mdot, vinf, rho_0, T, mu):
	return calcBubbleRadiusPRB(mdot, vinf, rho_0, calcConfineTimePRB(mdot, vinf, rho_0, T, mu))

###############################
### Adiabatic Bubble (FW)
###############################

def calcBubbleRadiusFAB(mdot, vinf, rho_0, time):
	L = calcL(mdot, vinf)
	return 0.88 * math.pow(L * (time**3) / rho_0, 0.2)

def calcBubbleRadiusPAB(mdot, vinf, rho_0, time):
	L = calcL(mdot, vinf)
	return 0.76 * math.pow(L * time * time * time / rho_0, 0.2)

def calcConfineTimeAB(mdot, vinf, rho_0, T, mu):
	L = calcL(mdot, vinf)
	cs = soundSpeed(T, mu)
	return 0.142 * math.sqrt(L / (rho_0 * (cs**5)))

def calcInjectionTimeAB(mdot, vinf, rho_0, R):
	L = calcL(mdot, vinf)
	return math.pow(R / (0.76 * math.pow(L / rho_0, 0.2)), 5.0 / 3.0)

def calcFinalRadiusFAB(mdot, vinf, rho_0, T, mu):
	return calcBubbleRadiusFAB(mdot, vinf, rho_0, calcConfineTimeAB(mdot, vinf, rho_0, T, mu))

def calcFinalRadiusPAB(mdot, vinf, rho_0, T, mu):
	return calcBubbleRadiusPAB(mdot, vinf, rho_0, calcConfineTimeAB(mdot, vinf, rho_0, T, mu))

def calcContactRadiusAB(mdot, vinf, rho_0, time, T, mu):
	L = calcL(mdot, vinf)
	cs = soundSpeed(T, mu)
	P0 = (cs**2) * rho_0
	return pow(3.0 * L * time / (10.0 * math.pi * P0), 1.0/3.0)

###############################
### Stromgren
###############################

def calcAlphaB(T):
	return 2.59e-13*math.pow(T / 10000.0, -0.7)

def calcStromgrenRadius(logQ, nh):
	alphaB = calcAlphaB(8000.0)
	return pow(3.0 * pow(10.0, logQ) / (4.0 * math.pi * nh * nh * alphaB), 1.0/3.0)

def calcSpitzerRadius(logQ, nh, t):
	RS = calcStromgrenRadius(logQ, nh)
	cs = soundSpeed(8000.0, 0.5)
	ts = RS / cs
	return RS * np.power(1.0 + np.sqrt(4.0 / 3.0) * ((7.0 * t)/(4.0 * ts)), 4.0 / 7.0)

def calcSpitzerRadius2(RS, t):
	cs = soundSpeed(8000.0, 0.5)
	ts = RS / cs
	return RS * np.power(1.0 + np.sqrt(4.0 / 3.0) * ((7.0 * t)/(4.0 * ts)), 4.0 / 7.0)

def calcStagnationRadius(logQ, nh):
	T = 8000.0
	ci = soundSpeed(8000.0, 0.5)
	co = soundSpeed(300.0, 1.0)
	cico = ci / co
	return math.pow(math.sqrt(4.0 / 3.0) * cico, 4.0 / 3.0) * calcStromgrenRadius(logQ, nh)

def calcStagnationTime(logQ, nh):
	ci = soundSpeed(8000.0, 0.5)
	co = soundSpeed(300.0, 1.0)
	cico = ci / co

	RS = calcStromgrenRadius(logQ, nh)
	ts = RS / ci

	return (4.0 / 7.0) * ((3.0 / 4.0)**0.5) * (((4.0 / 3.0)**(7.0 / 6.0)) * (cico**(7.0 / 3.0)) - 1.0) * ts

def calcStagnationRadius2(RS):
	T = 8000.0
	ci = soundSpeed(8000.0, 0.5)
	co = soundSpeed(300.0, 1.0)
	cico = ci / co
	return math.pow(math.sqrt(4.0 / 3.0) * cico, 4.0 / 3.0) * RS
