import argparse
import glob
import numpy as np
import scipy.interpolate
import linecache
import math

from torch_plot import TorchCFD, TorchInterpolatedGrid, TorchPlotter

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
	torchData = TorchCFD(inputfile, axial=True)
	t = torchData.t
	hii = torchData.get_var('hii')
	
	###    Interpolation set up.
	
	hiii = torchData.interpolate(hii, 'nearest')
	xi = torchData.xi[0]
	yi = torchData.xi[1]
	
	### Find ionisation radius.
	index = 0
	for i in range(0, nx):
		if hiii[i,0] < 0.5:
			index = i
			break

	if index != 0:
		fracLeft = hiii[index-1,0];
		interp = (0.5-hiii[index-1,0])/(hiii[index,0]-hiii[index-1,0]);
		IF = (xi[0, index-1] + interp)*(xi[0, index-1] + interp)
		if torchData.nd > 1:
			IF += yi[0, 0]*yi[0, 0]
	else:
		IF = 0;
		
	time.append(t)
	radius.append(math.sqrt(IF))

zipped = zip(time, radius)
zipped.sort()
time, radius = zip(*zipped)

out = open('if.csv', 'w')
for i in range(0, len(time)):
    out.write(str(time[i]) + ',' + str(radius[i]) + '\n')
out.close()
    
