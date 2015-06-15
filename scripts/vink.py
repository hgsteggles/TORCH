#Given a stellar mass, T_eff, logL and z, determine m_dot and v_inf
#using the Vink et al. (1999,2000,2001) prescriptions

import math
import sys
import argparse

print ("")

###	Parse arguements
parser = argparse.ArgumentParser(description='Calculates wind parameters.')
parser.add_argument('mass', metavar='mass', type=float, help='Star Mass.')
parser.add_argument('T_eff', metavar='T_eff', type=float, help='Effective Temperature.')
parser.add_argument('logL', metavar='logL', type=float, help='Logarithm of the Luminosity.')
parser.add_argument('radius', metavar='radius', type=float, help='Star Radius.')
parser.add_argument('-Z', metavar='Z', type=float, help='Metallicity.')

args = parser.parse_args()

mass = args.mass
T_eff = args.T_eff
logL = args.logL
radius = args.radius
Z = 0.02
if (args.Z != None):
	Z = args.Z
Z = Z/0.02

msol = 1.9891e33 # g
lsol = 3.846e33 # ergs.s-1
rsol = 6.955e10 # cm
G = 6.67384e-8 # cm3.g-1.s-2
stefb = 5.6704e-5 # erg.cm-2.s-1.K-4
yr = 3.15569e7 # s

logZ = math.log10(Z)
luminosity = math.pow(10.0, logL)
radius = radius*rsol
sigmae = 0.325
Gamma = 7.66*1.0e-5*sigmae*luminosity/mass
Gamma = min(Gamma, 0.999)
mass_eff = mass*(1.0 - Gamma)
#radius = math.sqrt(luminosity*lsol/(4.0*math.pi*stefb*T_eff*T_eff*T_eff*T_eff))
v_esc = math.sqrt(2.0*G*mass_eff*msol/radius)
v_esc *= 1.0e-5 # convert from cm.s-1 to km.s-1

rho = -14.94 + 3.1857*Gamma + 0.85*logZ # Eq. (23) from Vink et al. (2001) 

T_eff1 = 61200.0 + 2590.0*rho # Eq. (15) from Vink et al. (2001)
T_eff2 = 100000.0 + 6000.0*rho # Eq. (6)  from Vink et al. (2000)

def exit_program():
	print ("######### INPUTS ###########")
	print ("T_eff1    = " + str(T_eff1) + " K")
	print ("T_eff2    = " + str(T_eff2) + " K")
	print ("mass      = " + str(mass) + " Mo")
	print ("T_eff     = " + str(T_eff) + " K")
	print ("log(L/Lo) = " + str(logL))
	print ("radius = " + str(radius/rsol) + " Ro")
	print ("Z         = " + str(Z) + " Zo")
	print ("############################")
	sys.exit()

ratio = 0.0

if T_eff1 <= T_eff2:
	print ("These stellar parameters are unrealistic. No mass loss rate can be calculated.")
	exit_program()

if T_eff < T_eff1:
	if T_eff < T_eff2:
		ratio = 0.7
		logMdot = 2.210*math.log10(luminosity/(1.0e5)) + 1.07*math.log10(T_eff/20000.0) + 0.85*logZ
		logMdot -= 5.990 + 1.339*math.log10(mass/30.0) + 1.601*math.log10(ratio/2.0)
	else:
		# Vink et al. (2000) Eq. 13
		ratio = 1.3
		logMdot = 2.210*math.log10(luminosity/(1.0e5)) + 1.07*math.log10(T_eff/20000.0) + 0.85*logZ
		logMdot -= 6.688 + 1.339*math.log10(mass/30.0) + 1.601*math.log10(ratio/2.0)

if T_eff > T_eff1:
	# Vink et al. (2000) Eq. 12
	ratio = 2.6
	logMdot = 2.194*math.log10(luminosity/1.0e5) + 0.933*math.log10(T_eff/40000.0) + 0.85*logZ
	logMdot -= 6.697 + 1.313*math.log10(mass/30.0) + 1.226*math.log10(ratio/2.0) + 10.92*math.pow(math.log10(T_eff/40000.0), 2.0)

if T_eff == T_eff1 or T_eff == T_eff2:
	print ("The star is exactly on a jump. No mass loss rate can be calculated.")
	exit_program()
	
mdot = math.pow(10.0, logMdot)
mdot *= msol/yr # convert from msol.yr-1 to g.s-1
v_inf = ratio*v_esc
#v_inf *= 1.0e5 # convert from km.s-1 to cm.s-1

if math.isinf(v_inf) or math.isnan(v_inf):
	print ("v_inf is NaN or infinity.")
	exit_program()
else:
	print ("######### INPUTS ###########")
	print ("T_eff1    = " + str(T_eff1) + " K")
	print ("T_eff2    = " + str(T_eff2) + " K")
	print ("mass      = " + str(mass) + " Mo")
	print ("T_eff     = " + str(T_eff) + " K")
	print ("log(L/Lo) = " + str(logL))
	print ("radius = " + str(radius/rsol) + " Ro")
	print ("Z         = " + str(Z))
	print ("")
	print ("######### OUTPUTS ###########")
	print ("mdot      = " + str(mdot) + " g/s")
	print ("v_inf     = " + str(v_inf) + " km/s")
	print ("")





