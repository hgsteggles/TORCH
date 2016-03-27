import matplotlib.pyplot as plt
import math

PC2CM = 3.09e18
G = 6.6726e-8
hydrogenMass = 1.674e-24
specificGasConstant = 8.314462e7

T = 100
nHI = 80000
alpha = 1
rc = 0.01*PC2CM
rc2 = rc*rc
RS = 0.35*PC2CM
n0 = nHI*math.pow(1 + RS*RS/rc2, alpha)

RMIN = 0.00*PC2CM
RMAX = 0.40*PC2CM
N = 100
DR = RMAX/float(N)

def den(r):
	return n0*hydrogenMass*math.pow(1 + (r*r)/rc2, -alpha)

def mass(r, alpha):
	if alpha == 1:
		return 4*math.pi*n0*hydrogenMass*rc2*(r - rc*math.atan(r/rc))
	else:
		raise ValueError('invalid value of alpha: ' + str(alpha))

r = [0] * N
for i  in range(0, N):
	r[i] = (0.5 + i)*DR + RMIN;

clump_mass = mass(RS, alpha)
pot_pm = [0] * N
for i  in range(0, N):
	pot_pm[i] = -G*clump_mass/r[i]

pot_cl = [0] * N
g_coeff = -2.0*n0*hydrogenMass*specificGasConstant*T*alpha/rc2
for i  in range(0, N):
	pot_cl[i] = g_coeff*r[i]*math.pow(1 + (r[i]*r[i])/rc2, -alpha-1.0)/den(r[i])

mass_req = -(RS/G)*g_coeff*RS*math.pow(1 + (RS*RS)/rc2, -alpha-1.0)/den(RS)
pot_rm = [0] * N
for i  in range(0, N):
	pot_rm[i] = -G*mass_req/r[i]

plt.plot(r, pot_pm)
plt.plot(r, pot_cl)

plt.show()
