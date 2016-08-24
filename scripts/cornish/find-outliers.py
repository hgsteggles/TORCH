import numpy as np
import matplotlib.pyplot as plt

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

simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey.txt", skip_header=1)

for i in range(len(simulated_survey[:,0])):
	if abs((simulated_survey[i,8] - simulated_survey[i,9]) / simulated_survey[i,9]) > 0.5:
		print int(simulated_survey[i,0])
