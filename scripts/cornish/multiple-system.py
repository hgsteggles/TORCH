import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("cdat", "../torchpack/cornishdata.py")
import cdat
cornish_data = cdat.CornishData()

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData(useSet4=True)
masses = mp1_data.masses

star_data = np.genfromtxt("data/cornish/starpop-hm-w-a-s.txt", skip_header=1)

indexes = []

for i in range(1, len(star_data[:,0]) + 1):
	if not cornish_data.exclusion_set[i]:
		indexes.append(i)

raw_sizes = np.genfromtxt("data/cornish/raw-size.txt", skip_header=1)
ang_sizes = np.zeros(722, dtype=np.float64)
phy_sizes = np.zeros(722, dtype=np.float64)
for i in range(len(raw_sizes[:,0])):
	ang_sizes[raw_sizes[i,0]-1] = raw_sizes[i,1]
	phy_sizes[raw_sizes[i,0]-1] = raw_sizes[i,2]

raw_tot_fluxes = np.genfromtxt("data/cornish/raw-total-flux.txt", skip_header=1)
tot_fluxes = np.zeros(722, dtype=np.float64)
for i in range(len(raw_tot_fluxes[:,0])):
	tot_fluxes[raw_tot_fluxes[i,0]-1] = raw_tot_fluxes[i,1]

n = 0

f = open("data/cornish/excl-multiple.txt", 'w')

for i in range(len(indexes)):
	for j in range(i + 1, len(indexes)):
		ind1 = indexes[i]
		ind2 = indexes[j]
		ldiff = star_data[ind1 - 1,16] - star_data[ind2 - 1,16]
		bdiff = star_data[ind1 - 1,15] - star_data[ind2 - 1,15]
		r = math.sqrt(ldiff**2 + bdiff**2) * 60.0 * 60.0

		if r < 1.5:
			n += 1
			print str(ind1) + " " + str(ind2) + " " + str(ang_sizes[ind1 - 1]) + " " + str(ang_sizes[ind2 - 1]) + " " + str(tot_fluxes[ind1 - 1]) + " " + str(tot_fluxes[ind2 - 1])



print n

f.close()