import numpy as np

###	Data set up.
data = np.genfromtxt("refdata/zams-stripped.txt", skip_header=1)

mass = data[:,0] # [Msun]
temp = data[:,1]
logQ = data[:,6] # log[Q / phot.s-1]
mdot = data[:,8]/6.341958397e25 # [Msun.yr-1]
vinf = data[:,9]/100000.0 # [km.s-1]
logL = data[:,2] # log[L / Lsun]
rstar = data[:,3] # [Rsun]

den = np.array([0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4])

for i in range(len(den)):
	for j in range(len(mass)):
		id = "{:2d}".format(9*i + j + 1)
		m = "{:3d}".format(int(mass[j]))
		t = "{:1.2f}".format(temp[j])
		r = "{:4.2f}".format(rstar[j])
		l = "{:4.2f}".format(logL[j])
		q = "{:4.2f}".format(logQ[j])
		md = "{:8.2e}".format(mdot[j])
		v = "{:8.2e}".format(vinf[j])
		d = "{:8.2e}".format(den[i])
		print id + " & " + m + " & " + t + " & " + r + " & " + l + " & " + q + " & " + md + " & " + v + " & " + d + " \\\\"

