import numpy as np
import math
from astropy.io import fits
import argparse
import multiprocessing
import signal

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

def getPeak(filename):
	hdu_list = fits.open(filename)
	data = hdu_list[0].data

	peak = 0
	for i in range(len(data[0,0,:,:])):
		for j in range(len(data[0,0,i,:])):
			peak = max(peak, data[0,0,i,j])


	return peak

def getLine(index):
	try:
		ipad = cornish_data.star_id_str(index)
		dirname = cornish_data.dirname + "/star_" + ipad

		peak_0 = getPeak(dirname + "/radio_0/intensity_casa_ff.fits")
		peak_1 = getPeak(dirname + "/radio_1/intensity_casa_ff.fits")

		mass = star_data[index - 1,0]
		imass = 0
		for j in range(len(masses)):
			if masses[j] > mass:
				break
			imass = j
		interp = (mass - masses[imass])/(masses[imass + 1] - masses[imass])

		peak = peak_0 + interp * (peak_1 - peak_0)

		return str(index) + " " + str(peak * 1000.0)
	except:
		return ""

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

if __name__ == '__main__':
	pool = multiprocessing.Pool(8, init_worker) #use all available cores, otherwise specify the number you want as an argument
	try:
		lines = pool.map(getLine, range(1, nstars + 1))
		ofile = open(cornish_data.dirname + "/peak-fluxes.txt", 'w')
		ofile.write("star_id flux[mJy.beam-1]\n")

		for line in lines:
			if line != "":
				ofile.write(line + '\n')

		ofile.close()
	except KeyboardInterrupt:
		print "Caught KeyboardInterrupt, terminating workers"
		pool.terminate()
		pool.join()

