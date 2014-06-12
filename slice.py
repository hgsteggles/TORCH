#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import sys, getopt
import math
import linecache
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import matplotlib.colors as mcolors
import sod

inputfile_hll = 'hll.txt'
inputfile_hllc = 'hllc.txt'
inputfile_rot = 'rot.txt'
outputfile = ''
paramsfile = ''
figformat = 'png'
var_type = 'den'
y_slice = 49
DPI = 300

def plot_solver(ifile):
	#Get data
	data = np.genfromtxt(ifile,skip_header=0)

	if data.shape[1] == 7:
		x = data[:,0]
		y = data[:,1]
		if var_type == 'den':
			vid = 2
		if var_type == 'pre':
			vid = 3
		if var_type == 'hii':
			vid = 4
		x = np.zeros(1)
		var = np.zeros(1)
		appending = False
		for i in range(1, data.shape[0]):
			if data[i,1] == y_slice:
				if appending == False:
					x[0] = data[i,0]
					var[0] = data[i,vid]
					appending = True
				else:
					x = np.append(x, data[i,0])
					var = np.append(var, data[i,vid])

		plt.xlim([x.min()-0.1*(x.max()-x.min()), x.max()+0.1*(x.max()-x.min())]);
		plt.ylim([var.min()-0.1*(var.max()-var.min()), var.max()+0.1*(var.max()-var.min())]);
		plt.scatter(x[::1],var[::1], marker='o',color='blue',s=20)

	else:
		print "Can not handle %d columns" % data.shape[1]
		return

def main():
	fig=plt.figure()

	plot_solver(inputfile_hll)
	plot_solver(inputfile_hllc)
	plot_solver(inputfile_rot)

	x_an = np.array([0, 49.5, 49.5, 99])
	y_an = np.array([1.4, 1.4, 1.0, 1.0])
	plt.plot(x_an,y_an, color='black')

	extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
	fig.savefig(outputfile, format=figformat, dpi=DPI, bbox_inches=extent, pad_inches=0)
	print 'Plotted ' + outputfile

main()


