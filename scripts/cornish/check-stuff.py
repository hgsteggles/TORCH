import numpy as np
import os.path

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

indexes = np.fromfile(cornish_data.dirname + "/valid-pre-casa.txt", dtype=int, sep='\n')

for index in indexes:
	ipad = cornish_data.star_id_str(index)

	for i in range(2):
		filename = cornish_data.dirname + "/star_" + ipad + "/radio_" + str(i) + "/intensity_casa_ff.fits"
		if not os.path.isfile(filename):
			print "usr/cornish_" + str(args.iden) + "/cornish_star_" + ipad + "_" + str(i) + ".py"
