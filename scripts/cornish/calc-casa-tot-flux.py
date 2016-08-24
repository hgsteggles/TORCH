import numpy as np
import math
from astropy.io import fits
import argparse
import re

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

pattern = re.compile("Integrated:   (.*)")

def getTotalFlux(filename):
	flux_string = ""
	for i, line in enumerate(open(filename)):
		for match in re.finditer(pattern, line):
			flux_string = match.groups()[0]

	scale = 1.0
	if flux_string.split()[3] == 'Jy':
		scale = 1000.0
	elif flux_string.split()[3] == 'uJy':
		scale = 1.0 / 1000.0

	return float(flux_string.split()[0]) * scale

remaining_casa = np.genfromtxt(cornish_data.dirname + "/valid-post-casa.txt", dtype=int)

ofile = open(cornish_data.dirname + "/casa-tot-fluxes.txt", 'w')
ofile.write("star_id ang_diameter [\"]\n")

for index in remaining_casa:
	ipad = cornish_data.star_id_str(index)
	dirname = cornish_data.dirname + "/star_" + ipad + "/"

	try:
		total_flux_0 = getTotalFlux(dirname+"radio_0/gaussfit-full.txt")
		total_flux_1 = getTotalFlux(dirname+"radio_1/gaussfit-full.txt")
	except:
		continue

	mass = star_data[index - 1,0]
	imass = 0
	for j in range(len(masses)):
		if masses[j] > mass:
			break
		imass = j
	interp = (mass - masses[imass])/(masses[imass + 1] - masses[imass])

	total_flux = total_flux_0 + interp * (total_flux_1 - total_flux_0)
	ofile.write(str(index) + " " + str(total_flux) + "\n")

ofile.close()
