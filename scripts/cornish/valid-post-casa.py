import numpy as np
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

remaining_set = np.zeros(len(star_data[:,0]) + 1, dtype=bool)

valid_pre_casa = np.genfromtxt(cornish_data.dirname + "/valid-pre-casa.txt", skip_header=0, dtype=int)
for id in valid_pre_casa:
	remaining_set[id] = True

weak_flux_ids = np.genfromtxt(cornish_data.dirname + "/excl-casa-weak-fluxes.txt", skip_header=1, dtype=int)
for id in weak_flux_ids:
	remaining_set[id] = False


ofile = open(cornish_data.dirname + "/valid-post-casa.txt", 'w')

for i in range(1, len(remaining_set)):
	if (remaining_set[i]):
		ofile.write(str(i) + '\n')

ofile.close()