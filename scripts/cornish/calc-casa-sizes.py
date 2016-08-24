import numpy as np
import math
from astropy.io import fits
import argparse

###	Parse arguements
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData(useSet4=True)
masses = mp1_data.masses

load_src("cdat", "../torchpack/cornishdata.py")
import cdat
cornish_data = cdat.CornishData(args.iden)
star_data = cornish_data.star_data
nstars = len(star_data[:,0])

def getAngSize(filename):
	gauss_fit = np.genfromtxt(filename, delimiter=', ', dtype=str)
	avg = (float(gauss_fit[3].split()[0]) + float(gauss_fit[4].split()[0])) / 2.0
	return math.sqrt(avg**2 - (1.5)**2)

remaining_casa = np.genfromtxt(cornish_data.dirname + "/valid-post-casa.txt", dtype=int)

ofile = open(cornish_data.dirname + "/casa-sizes.txt", 'w')
ofile.write("star_id ang_diameter [\"]\n")

for index in remaining_casa:
	ipad = cornish_data.star_id_str(index)
	dirname = cornish_data.dirname + "/star_" + ipad + "/"

	try:
		angsize_0 = getAngSize(dirname+"radio_0/gaussfit.txt")
		angsize_1 = getAngSize(dirname+"radio_1/gaussfit.txt")
	except:
		continue

	mass = star_data[index - 1,0]
	imass = 0
	for j in range(len(masses)):
		if masses[j] > mass:
			break
		imass = j
	interp = (mass - masses[imass])/(masses[imass + 1] - masses[imass])

	angsize = angsize_0 + interp * (angsize_1 - angsize_0)
	ofile.write(str(index) + " " + str(angsize) + "\n")

ofile.close()
