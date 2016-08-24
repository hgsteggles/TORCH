import sys
import warnings
import numpy as np

###	Parse arguements
import argparse
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

load_src("cdat", "../torchpack/cornishdata.py")
import cdat

nhs = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]
masses = [6.0, 9.0, 12.0, 15.0, 20.0, 30.0, 40.0, 70.0, 120.0]

ofile = open("data/latitude_" + str(args.iden) + "/final-survey.txt", 'w')
ofile.write("star_id mcur[Msun] t_ms[kyr] age[kyr] d_sun[kpc] phy_size[pc] ang_size[\"] casa_ang_size[\"] tot_flux[mJy] casa_tot_flux[mJy] casa_pk_flux[mJy.beam-1] g_long[deg] g_lat[deg]\n")

cornish_data = cdat.CornishData(args.iden)
star_data = cornish_data.star_data
nstars = len(star_data[:,0])

### Data
content = open(cornish_data.dirname + "/final-survey" + ".txt").readlines()
survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + ".txt", skip_header=1)
cornish_survey = np.genfromtxt("data/cornish/cornish-uchiis.txt", skip_header=1, delimiter=',')

xbins = np.arange(10, 65, 2)
counts, xedges = np.histogram(cornish_survey[:,1], bins=xbins)

sim_counts = np.zeros(counts.shape, dtype=float)

for i in range(len(survey[:,11])):
	l = survey[i,11]

	for j in range(len(counts)):
		if l >= xedges[j] and l < xedges[j + 1]:
			if sim_counts[j] < counts[j]:
				sim_counts[j] = sim_counts[j] + 1
				ofile.write(content[i + 1])
			break

ofile.close()
