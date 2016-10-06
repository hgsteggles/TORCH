import argparse
import glob
import numpy as np
import scipy.interpolate
import linecache
import math

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
import torch

###    Parse arguements
parser = argparse.ArgumentParser(description='Plots wind radius vs time.')
parser.add_argument('input_dir', metavar='input_directory', type=str,
				help='Input directory with CFD output data.')

args = parser.parse_args()

input_dir = args.input_dir

if not input_dir.endswith('/'):
	input_dir += '/'
input_regex = input_dir + '*.txt'

time = []
radius = []

for inputfile in glob.iglob(input_regex):
	### Open file.
	f = open(inputfile,'r')
	### Get Data
	torchData = torch.CFD_Data(inputfile, axial=False)
	t = torchData.t
	hii = torchData.get_var_raw('hii')
	
	###    Interpolation set up.

	xi = torchData.xi[0]
	
	### Find ionisation radius.
	index = 0
	for i in range(0, len(hii)):
		if hii[i] < 0.5:
			index = i
			break

	if index != 0:
		fracLeft = hii[index - 1]
		interp = (0.5 - hii[index - 1]) / (hii[index] - hii[index - 1])
		IF = (xi[index - 1] + interp * (xi[index] - xi[index - 1]))
	else:
		interp = 0
		IF = 0

	print inputfile + " " + str(index) + " " + str(xi[index] * 3.09e18) + " " + str(interp)

	time.append(t)
	radius.append(IF)

zipped = zip(time, radius)
zipped.sort()
time, radius = zip(*zipped)

out = open('if.csv', 'w')
for i in range(0, len(time)):
    out.write(str(time[i]) + ',' + str(radius[i]) + '\n')
out.close()

