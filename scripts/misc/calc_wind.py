import argparse
import glob
import numpy as np
import scipy.interpolate
import linecache
from torch_plot import TorchCFD

###	Parse arguements
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
	### Data set up.
	torchData = TorchCFD(inputfile, axial=True)
	den = torchData.get_var('den')
	deni = torchData.interpolate(den, 'nearest')

	### Find bubble radius.
	index = 0
	max_den = 0
	for i in range(11, torchData.nx):
		if deni[i,0] > max_den:
			index = i
			max_den = deni[i,0]
	time.append(t)
	radius.append(torchData.xi[0, index])

zipped = zip(time, radius)
zipped.sort()
time, radius = zip(*zipped)

out = open('wind.csv', 'w')
for i in range(0, len(time)):
	out.write(str(time[i]) + ',' + str(radius[i]) + '\n')
out.close()
	
