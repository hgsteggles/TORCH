import numpy as np

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

data = np.genfromtxt(cornish_data.dirname + "/sizes.txt", skip_header=1)

remaining_set = np.zeros(nstars + 1, dtype=bool)

for i in range(0, len(data[:,0])):
	index = int(data[i,0])
	ang_diameter = data[i,1]
	phy_diameter = data[i,2]

	if ang_diameter <= 24 and phy_diameter <= 1.0:
		remaining_set[index] = True


ofile = open(cornish_data.dirname + "/excl-large-sizes.txt", 'w')
ofile.write("star_id [(ang_size > 24) || (phy_diameter > 1.0 pc)]\n")

for index in range(1, nstars + 1):
	if not remaining_set[index]:
		ofile.write(str(int(index)) + "\n")

ofile.close()
