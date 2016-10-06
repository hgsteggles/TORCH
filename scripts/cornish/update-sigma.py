import numpy as np
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

sigma = 7

def ga2equ(ga):
    l = math.radians(ga[0])
    b = math.radians(ga[1])

    # North galactic pole (J2000) -- according to Wikipedia
    pole_ra = math.radians(192.859508)
    pole_dec = math.radians(27.128336)
    posangle = math.radians(122.932-90.0)

    # North galactic pole (B1950)
    #pole_ra = radians(192.25)
    #pole_dec = radians(27.4)
    #posangle = radians(123.0-90.0)

    ra = math.atan2( (math.cos(b)*math.cos(l-posangle)), (math.sin(b)*math.cos(pole_dec) - math.cos(b)*math.sin(pole_dec)*math.sin(l-posangle)) ) + pole_ra
    dec = math.asin( math.cos(b)*math.cos(pole_dec)*math.sin(l-posangle) + math.sin(b)*math.sin(pole_dec) )

    return np.array([math.degrees(ra), math.degrees(dec)])

data = np.genfromtxt(cornish_data.dirname + "/final-survey-raw.txt", skip_header=1)

with open(cornish_data.dirname + "/final-survey-raw.txt") as f:
	content = f.readlines()

ofile = open(cornish_data.dirname + "/final-survey.txt", 'w')
ofile.write(content[0])

for i in range(len(data[:,0])):
	casa_peak_flux = data[i,10]
	long = data[i,11]
	lat = data[i,12]

	ra_dec = ga2equ([lat, long])

	if (ra_dec[1] >= 14.2 and casa_peak_flux > 0.25 * sigma) or (ra_dec[1] < 14.2 and casa_peak_flux > 0.35 * sigma):
		ofile.write(content[i + 1])

ofile.close()