import sys
import warnings
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

### Data
simulated_survey = np.genfromtxt(cornish_data.dirname + "/final-survey.txt", skip_header=1)
cornish_survey = np.genfromtxt("data/cornish/cornish-uchiis.txt", skip_header=1, delimiter=',')

print "cornish mean flux   = " + str(sum(cornish_survey[:,9]) / float(len(cornish_survey[:,9]))) + " mJy"
print "simulated mean flux = " + str(sum(simulated_survey[:,9]) / float(len(simulated_survey[:,9]))) + " mJy"
