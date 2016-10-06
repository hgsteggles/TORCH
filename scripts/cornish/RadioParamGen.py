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

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData(useSet4=True)
masses = mp1_data.masses

load_src("cdat", "../torchpack/cornishdata.py")
import cdat
cornish_data = cdat.CornishData(args.iden)
cornish_data.dirname = "data/cornish_examples"
star_data = cornish_data.star_data

nstars = len(star_data[:,0])

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

def write_file(ist, im, id, it, isLower, l, b, dist, theta):
	data_dirname = mp1_data.getDataDirname(im, id, it)
	data_filename = mp1_data.getDataFilename(im, id, it)
	index = mp1_data.getIndex(im, id)
	config_id = 0 if isLower else 1
	istpad = cornish_data.star_id_str(ist)
	ra_dec = ga2equ([l, b])

	f = open(cornish_data.dirname + "/star_" + istpad + "/radioconfig_" + str(config_id) + ".lua", 'w')

	string = """-- RadioRT Parameters
Parameters = {
	torch_params_filename =      "../../Torch/""" + data_dirname  + """/parameters_""" + str(index) + """.lua",
	torch_data_filename =        "../../Torch/""" + data_filename + """",
	output_directory =           "../../Torch/""" + cornish_data.dirname + """/star_""" + istpad + """/radio_""" + str(config_id) + """/",

	sampling =                   4.0,
	dopplerShifted =             false,
	dopp_shift_phi_inc =         0.01,

	rbeam_degrees =              (3.4/2.0)/3600.0,
	distance =                   """ + str(dist) + """,
	right_ascension =            """ + str(ra_dec[0]) + """,
	declination =                """ + str(ra_dec[1]) + """,
	theta =                      """ + str(theta) + """,
	phi =                        0.0,

	frequency =                  5.0e9,
	bandwidth =                  0,
	nchannels =                  1,
	nlevel =                     0,
	stokes =                     1,

	turb_broadening =            1.0e6,
	vLOS =                       0,

	integratingFF =              true,
	integratingRL =              false,
	resolution_scale =           2,
}"""
	f.write(string)

for i in range(nstars):
	l = star_data[i,15]
	b = star_data[i,16]
	dist = star_data[i,4]

	mass = star_data[i,0]
	imass = 0
	for j in range(len(masses)):
		if masses[j] > mass:
			break
		imass = j

	age = star_data[i,2]

	if (int)(age / 2.0 + 0.5) > 100:
		if args.iden < 3:
			continue

	itime = min(100, (int)(age / 2.0 + 0.5))

	theta = random.uniform(-180.0, 180.0)
	if theta > 90.0:
		theta = 180.0 - theta
	if theta < -90.0:
		theta = -180.0 - theta

	write_file(i + 1, imass, args.iden, itime, True, l, b, dist, theta)
	write_file(i + 1, imass + 1, args.iden, itime, False, l, b, dist, theta)

	print mp1_data.getDataFilename(imass, args.iden, itime) + ".gz"
	print mp1_data.getDataFilename(imass + 1, args.iden, itime) + ".gz"

	index_0 = mp1_data.getIndex(imass, args.iden)
	index_1 = mp1_data.getIndex(imass + 1, args.iden)

	padi = cornish_data.star_id_str(i + 1)

	#if index_0 == 44:
	#	print "../Torch/data/cornish_" + str(args.iden) + "/star_" + ipad + "/radioconfig_0.lua"
	#elif index_1 == 44:
	#	print "../Torch/data/cornish_" + str(args.iden) + "/star_" + ipad + "/radioconfig_1.lua"

	#if not mp1_data.isSet2(imass, iden, itime):
	#	print str(i + 1) + " in set4"

	#if imass == 4:
	#	if not mp1_data.isSet2(imass, iden, itime):
	#		print str(i) + " uses data_23"
	#elif imass == 5:
	#	if not mp1_data.isSet2(imass, iden, itime):
	#		print str(i) + " uses data_24"
