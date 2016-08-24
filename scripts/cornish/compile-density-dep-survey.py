import sys
import warnings
import numpy as np

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

load_src("cdat", "../torchpack/cornishdata.py")
import cdat

surveys = []
nhs = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]
masses = [6.0, 9.0, 12.0, 15.0, 20.0, 30.0, 40.0, 70.0, 120.0]

alpha = 5.0

a = (nhs[-1] - nhs[0]) / (masses[-1]**alpha - masses[0]**alpha)
b = nhs[0] - a * masses[0]**alpha

def nh_from_mass(mass):
	return a * mass**alpha + b

def index_from_nh(nh):
	hi = 0
	lo = len(nhs) - 1

	for i in range(len(nhs) - 1, -1, -1):
		if nh < nhs[i]:
			hi = i

	for i in range(len(nhs)):
		if nh >= nhs[i]:
			lo = i

	if lo == hi or hi == 0:
		return hi
	elif lo == len(nhs) - 1:
		return lo
	elif abs(nhs[lo] - nh) < abs(nhs[hi] - nh):
		return lo
	else:
		return hi

ofile = open("data/cornish-ne2001/cornish-tms/density_4/final-survey.txt", 'w')
ofile.write("star_id mcur[Msun] t_ms[kyr] age[kyr] d_sun[kpc] phy_size[pc] ang_size[\"] casa_ang_size[\"] tot_flux[mJy] casa_tot_flux[mJy] casa_pk_flux[mJy.beam-1] g_long[deg] g_lat[deg]\n")

for i in range(5):
	cornish_data = cdat.CornishData(i)
	star_data = cornish_data.star_data
	nstars = len(star_data[:,0])

	### Data
	content = open(cornish_data.dirname + "/final-survey" + ".txt").readlines()
	survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + ".txt", skip_header=1)

	for j in range(len(survey[:,0])):
		mcur = survey[j,1]
		if i == index_from_nh(nh_from_mass(mcur)):
			ofile.write(content[j + 1])

ofile.close()
