import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import math
import multiprocessing
import signal

###	Parse arguements
import argparse
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("cdat", "../torchpack/cornishdata.py")
import cdat
cornish_data = cdat.CornishData(args.iden)

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData(useSet4=True)
masses = mp1_data.masses

star_data = cornish_data.star_data
nstars = len(star_data[:,0])

def nPixels(data, mpix):
	n = 0
	for i in range(len(data)):
		for j in range(len(data[i])):
			if data[i][j] > mpix / 1.0e11:
				n += 1
	return n

def getAngSizeArcSec(filename):
	hdu_list = fits.open(filename)
	pixdeg = hdu_list[0].header['CDELT2']
	pixas = pixdeg * 60.0 * 60.0
	mpix = hdu_list[0].header['MPIX']
	npixels = nPixels(hdu_list[0].data, mpix)

	return 2.0 * np.sqrt(npixels * pixas * pixas / math.pi)

def getPhysSizeParsec(filename):
	hdu_list = fits.open(filename)
	pixdeg = hdu_list[0].header['CDELT2']
	pixrad = pixdeg * math.pi / 180.0
	dist = hdu_list[0].header['DIST']
	pixpc = pixrad * dist * 1000.0
	mpix = hdu_list[0].header['MPIX']
	npixels = nPixels(hdu_list[0].data, mpix)

	return 2.0 * np.sqrt(npixels * pixpc * pixpc / math.pi)

def getLine(index):
	ipad = cornish_data.star_id_str(index)
	dirname = cornish_data.dirname + "/star_" + ipad

	try:
		angsize_0 = getAngSizeArcSec(dirname + "/radio_0/intensity_pixel_ff.fits")
		angsize_1 = getAngSizeArcSec(dirname + "/radio_1/intensity_pixel_ff.fits")
		physize_0 = getPhysSizeParsec(dirname + "/radio_0/intensity_pixel_ff.fits")
		physize_1 = getPhysSizeParsec(dirname + "/radio_1/intensity_pixel_ff.fits")
	except:
		return ""

	mass = star_data[index - 1,0]
	imass = 0
	for j in range(len(masses)):
		if masses[j] > mass:
			break
		imass = j
	interp = (mass - masses[imass])/(masses[imass + 1] - masses[imass])

	angsize = angsize_0 + interp * (angsize_1 - angsize_0)
	physize = physize_0 + interp * (physize_1 - physize_0)

	line = str(index) + " " + str(angsize) + " " + str(physize)

	return line

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

if __name__ == '__main__':
	pool = multiprocessing.Pool(8, init_worker) #use all available cores, otherwise specify the number you want as an argument
	try:
		lines = pool.map(getLine, range(1, nstars + 1))
	except KeyboardInterrupt:
		print "Caught KeyboardInterrupt, terminating workers"
		pool.terminate()
		pool.join()

ofile = open(cornish_data.dirname + "/sizes.txt", 'w')
ofile.write("star_id ang_diameter[\"] phy_diameter[pc]\n")

for line in lines:
	if line != "":
		ofile.write(line + '\n')

ofile.close()
