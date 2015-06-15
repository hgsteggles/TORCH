import numpy as np
import matplotlib.pyplot as plt


data = np.genfromtxt("hummer.txt",skip_header=1)
logT = data[:,0]
T = np.power(10.0, logT)
aB = data[:,2]/np.sqrt(T)
aB2 = 2.59e-13*np.power(T/10000.0, -0.7);

print data[:,2]

aB = np.log10(aB)
aB2 = np.log10(aB2)

plt.plot(T, aB, label='hummer')
plt.plot(T, aB2, label='analyt')
plt.xlim([0, 10000])
plt.legend()

plt.show()
