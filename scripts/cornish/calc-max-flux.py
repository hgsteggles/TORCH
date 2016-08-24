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

def getMaxFlux(filename):
	hdu_list = fits.open(filename)
	return hdu_list[0].header['MPIX']

ofile = open(cornish_data.dirname + "/max-fluxes.txt", 'w')
ofile.write("star_id max_flux [mJy.beam-1]\n")

for i in range(1, len(star_data[:,0]) + 1):
	ipad = "%03d" % (i,)
	dirname = "data/cornish/star_" + ipad + "/"

	max_flux_0 = getMaxFlux(dirname+"radio_0/intensity_beam_ff.fits")
	max_flux_1 = getMaxFlux(dirname+"radio_1/intensity_beam_ff.fits")

	mass = star_data[i - 1,0]
	imass = 0
	for j in range(len(masses)):
		if masses[j] > mass:
			break
		imass = j
	interp = (mass - masses[imass])/(masses[imass + 1] - masses[imass])

	max_flux = max_flux_0 + interp * (max_flux_1 - max_flux_0)
	ofile.write(str(i) + " " + str(max_flux) + "\n")

ofile.close()