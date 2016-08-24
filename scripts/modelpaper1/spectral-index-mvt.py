import numpy as np
import math
from astropy.io import fits

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData()

print_sizes = False
print_mvd_table = False
print_mvt_table = True

data = np.genfromtxt("data/model_paper1/refdata/zams.txt", skip_header=1)

massList = data[:,0]
logQList = data[:,1]
mdotList = data[:,2]
vinfList = data[:,3]

CM2PC = 3.24078e-19
S2YR = 3.17098e-8
YR2S = 3.154e7
mh = 1.6738232e-24
R_GAS = 8.3144598e7

snapshots = [10, 20, 30, 40, 50]

def getSpectralIndices():
	sindex = np.zeros((5, 9))

	for col in range(5):
		for row in range(9):

			hdu_list0 = fits.open(mp1_data.getRadioDirname(row, 2, snapshots[col], 45, 1.4) + "/intensity_beam_ff.fits")
			hdu_list1 = fits.open(mp1_data.getRadioDirname(row, 2, snapshots[col], 45, 5) + "/intensity_beam_ff.fits")
			data0 = hdu_list0[0].data
			data1 = hdu_list1[0].data

			#data0[data0 < 0.4] = 0.4
			#data1[data1 < 0.4] = 0.4

			sindex[col][row] = (np.log10(np.sum(data1) / 1000.0) - np.log10(np.sum(data0) / 1000.0)) / ((np.log10(5.0) - np.log10(1.4)))

	return sindex


def printTable(sindices):
	for row in range(9):
		tableline = str(massList[row])
		for col in range(5):
			tableline = tableline + " & "
			if sindices[col][row] != None and not np.isnan(sindices[col][row]):
				tableline = tableline + str(sindices[col][row])
			else:
				tableline = tableline + "{\\textemdash}"

		print tableline + " \\\\"

		if row != 8:
			print "& & & & & \\\\"


printTable(getSpectralIndices())
















