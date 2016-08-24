import os.path

import numpy as np
import math
import random

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

def check_file(ist, isLower):
	config_id = 0 if isLower else 1
	istpad = "%03d" % (ist,)

	if not os.path.isfile(cornish_data.dirname + "/star_" + istpad + "/radio_" + str(config_id) + "/intensity_pixel_ff.fits"):
		print "../../Torch/" + cornish_data.dirname + "/star_" + istpad + "/radioconfig_" + str(config_id) + ".lua"



for i in range(nstars):
	check_file(i + 1, True)
	check_file(i + 1, False)
