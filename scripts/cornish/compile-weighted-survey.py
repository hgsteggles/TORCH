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

cornish_data = cdat.CornishData(0)
star_data = cornish_data.star_data
nstars = len(star_data[:,0])

weight_id = 2
alpha = -1.0

A = 0.0
for nh in nhs:
	A += nh**alpha
A = 1.0 / A

ofile = open("data/weights_" + str(weight_id) + "/weights.txt", 'w')
ofile.write("alpha = " + str(alpha) + "\n\n")

weights = []
for i in range(len(nhs)):
	weights.append(A * (nhs[i]**alpha))
	ofile.write(str(weights[i]) + "\n")

#weights = [0.3, 0.25, 0.2, 0.15, 0.1]

ofile.close()

cum_weights = [weights[0]]
for i in range(1, len(weights)):
	cum_weights.append(cum_weights[i-1] + weights[i])

bounds = []
for i in range(len(cum_weights)):
	if i == 0:
		bounds.append([0, (int)(cum_weights[i] * nstars)])
	else:
		bounds.append([(int)(cum_weights[i - 1] * nstars), (int)(cum_weights[i] * nstars)])

ofile = open("data/weights_" + str(weight_id) + "/final-survey.txt", 'w')
ofile.write("star_id mcur[Msun] t_ms[kyr] age[kyr] d_sun[kpc] phy_size[pc] ang_size[\"] casa_ang_size[\"] tot_flux[mJy] casa_tot_flux[mJy] casa_pk_flux[mJy.beam-1] g_long[deg] g_lat[deg]\n")

for i in range(5):
	cornish_data = cdat.CornishData(i)

	### Data
	content = open(cornish_data.dirname + "/final-survey" + ".txt").readlines()
	survey = np.genfromtxt(cornish_data.dirname + "/final-survey" + ".txt", skip_header=1)

	for j in range(len(survey[:, 0])):
		istar = survey[j, 0]
		if istar > bounds[i][0] and istar <= bounds[i][1]:
			ofile.write(content[j + 1])

ofile.close()
