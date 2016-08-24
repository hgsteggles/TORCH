import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate


data = np.genfromtxt("refdata/hummer.txt",skip_header=1)
logT = data[:,0]
T = np.power(10.0, logT)

bB = data[:,4] / np.sqrt(T)
bB_Spline = scipy.interpolate.interp1d(T, bB, kind='cubic')

Tnew = np.linspace(0, 1.0e8, num=41, endpoint=True)

bB_Cool = Tnew * bB_Spline(Tnew)

plt.plot(Tnew, bB_Cool, )

aB = data[:,2]/np.sqrt(T)
aB2 = 2.59e-13*np.power(T/10000.0, -0.7)

print data[:,2]

aB = np.log10(aB)
aB2 = np.log10(aB2)

plt.plot(T, aB, label='hummer')
plt.plot(T, aB2, label='analyt')
plt.xlim([0, 1.0e7])
plt.legend()

plt.show()
