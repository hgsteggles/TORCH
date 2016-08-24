import math
from astropy.io import fits

###	Parse arguements
import argparse
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

def getTotalFlux(filename):
	hdu_list = fits.open(filename)
	return hdu_list[0].header['TOTAL']

ofile = open(cornish_data.dirname + "/total-fluxes.txt", 'w')
ofile.write("star_id total_flux [mJy]\n")

for i in range(1, len(star_data[:,0]) + 1):
	ipad = cornish_data.star_id_str(i)
	dirname = cornish_data.dirname + "/star_" + ipad + "/"

	try:
		tot_flux_0 = getTotalFlux(dirname+"radio_0/intensity_pixel_ff.fits")
		tot_flux_1 = getTotalFlux(dirname+"radio_1/intensity_pixel_ff.fits")
	except:
		continue

	mass = star_data[i - 1,0]
	imass = 0
	for j in range(len(masses)):
		if masses[j] > mass:
			break
		imass = j
	interp = (mass - masses[imass])/(masses[imass + 1] - masses[imass])

	tot_flux = tot_flux_0 + interp * (tot_flux_1 - tot_flux_0)
	ofile.write(str(i) + " " + str(tot_flux) + "\n")

ofile.close()
