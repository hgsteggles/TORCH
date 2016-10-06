import numpy as np
import sys
import matplotlib.pyplot as plt
from astropy.io import fits
import math

# Parse arguments.
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
star_data = cornish_data.star_data
nstars = len(star_data[:,0])

remaining_set = np.zeros(nstars + 1, dtype=bool)

remaining_casa = np.genfromtxt(cornish_data.dirname + "/casa-sizes.txt", dtype=int, skip_header=1)[:,0]
for index in remaining_casa:
	remaining_set[index] = True

mcur = star_data[:,0]
t_ms = star_data[:,2]
age = star_data[:,3]
d_sun = star_data[:,4]
phy_size = np.zeros(nstars, dtype=float)
ang_size = np.zeros(nstars, dtype=float)
casa_ang_size = np.zeros(nstars, dtype=float)
tot_flux = np.zeros(nstars, dtype=float)
casa_tot_flux = np.zeros(nstars, dtype=float)
casa_peak_flux = np.zeros(nstars, dtype=float)
g_long = star_data[:,15]
g_lat = star_data[:,16]

sizes = np.genfromtxt(cornish_data.dirname + "/sizes.txt", skip_header=1)
for i in range(len(sizes)):
	index = int(sizes[i,0])
	ang_size[index - 1] = sizes[i,1]
	phy_size[index - 1] = sizes[i,2]

casa_sizes = np.genfromtxt(cornish_data.dirname + "/casa-sizes.txt", skip_header=1)
for i in range(len(casa_sizes)):
	index = int(casa_sizes[i,0])
	casa_ang_size[index - 1] = casa_sizes[i,1]

peak_fluxes = np.genfromtxt(cornish_data.dirname + "/peak-fluxes.txt", skip_header=1)
for i in range(len(peak_fluxes)):
	index = int(peak_fluxes[i,0])
	casa_peak_flux[index - 1] = peak_fluxes[i,1]

total_fluxes = np.genfromtxt(cornish_data.dirname + "/total-fluxes.txt", skip_header=1)
for i in range(len(total_fluxes)):
	index = int(total_fluxes[i,0])
	tot_flux[index - 1] = total_fluxes[i,1]

casa_total_fluxes = np.genfromtxt(cornish_data.dirname + "/casa-tot-fluxes.txt", skip_header=1)
for i in range(len(casa_total_fluxes)):
	index = int(casa_total_fluxes[i,0])
	casa_tot_flux[index - 1] = casa_total_fluxes[i,1]

ofile1 = open(cornish_data.dirname + "/final-survey-raw.txt", 'w')
ofile2 = open(cornish_data.dirname + "/final-survey.txt", 'w')
ofile1.write("star_id mcur[Msun] t_ms[kyr] age[kyr] d_sun[kpc] phy_size[pc] ang_size[\"] casa_ang_size[\"] tot_flux[mJy] casa_tot_flux[mJy] casa_pk_flux[mJy.beam-1] g_long[deg] g_lat[deg]\n")
ofile2.write("star_id mcur[Msun] t_ms[kyr] age[kyr] d_sun[kpc] phy_size[pc] ang_size[\"] casa_ang_size[\"] tot_flux[mJy] casa_tot_flux[mJy] casa_pk_flux[mJy.beam-1] g_long[deg] g_lat[deg]\n")

for i in range(1, nstars + 1):
	if (remaining_set[i]):
		print_str = str(i) + " " + str(mcur[i-1]) + " " + str(t_ms[i-1]) +  " " + str(age[i-1]) \
					+ " " + str(d_sun[i-1])\
					+ " " + str(phy_size[i-1]) + " " + str(ang_size[i-1])\
					+ " " + str(casa_ang_size[i-1]) + " " + str(tot_flux[i-1])\
					+ " " + str(casa_tot_flux[i-1]) + " " + str(casa_peak_flux[i-1])\
					+ " " + str(g_long[i-1]) + " " + str(g_lat[i-1]) + '\n'
		ofile1.write(print_str)
		ofile2.write(print_str)

ofile1.close()
ofile2.close()
