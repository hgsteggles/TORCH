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

remaining_set = np.ones(len(star_data[:,0]) + 1, dtype=bool)

weak_flux_ids = np.genfromtxt(cornish_data.dirname + "/excl-weak-fluxes.txt", skip_header=1, dtype=int)

for id in weak_flux_ids:
	remaining_set[id] = False

small_size_ids = np.genfromtxt(cornish_data.dirname + "/excl-large-sizes.txt", skip_header=1, dtype=int)
for id in small_size_ids:
	remaining_set[id] = False

ofile = open(cornish_data.dirname + "/valid-pre-casa.txt", 'w')

for i in range(1, len(remaining_set)):
	if (remaining_set[i]):
		ofile.write(str(i) + '\n')

ofile.close()