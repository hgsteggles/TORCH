import sys
import warnings
import numpy as np
from scipy import stats

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('dir', metavar='dir', type=str, help='Directory containting galaxy population data.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

### Data
simulated_survey = np.genfromtxt(args.dir + "/starpopfinal.txt", skip_header=1)

mcur = simulated_survey[:, 0]
mfin = simulated_survey[:, 1]

m = mfin

xbins = np.logspace(np.log10(6.0), np.log10(70), 40)
#xbins = np.arange(6.0, 70.0, 20)

counts_0, xedges_0 = np.histogram(m, bins=xbins)

x = []
for i in range(len(xedges_0) - 1):
	x.append(np.sqrt(xedges_0[i] * xedges_0[i + 1]))

for i in range(len(counts_0)):
	counts_0[i] = counts_0[i] / (xedges_0[i + 1] - xedges_0[i])

slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(x), np.log10(counts_0))

print "slope = " + str(slope)