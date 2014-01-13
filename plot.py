#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import sys, getopt
import math
import linecache
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext

import sod

inputfile = ''
outputfile = ''
paramsfile = ''
figformat = 'png'
zlog = True

def parseArg(argv):
	global inputfile
	global outputfile
	try:
		opts, args = getopt.getopt(argv,"hi:")
	except getopt.GetoptError:
		print 'ERROR: plot.py -i <inputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'plot.py -i <inputfile> '
			sys.exit()
		elif opt in ("-i"):
			inputfile = arg
	if inputfile == '':
		print 'plot.py -i <inputfile>'
		sys.exit(2)
	elif not '.' in inputfile:
		inputfile += '.xq'
	outputfile = inputfile.partition('.')[0] + '.' + figformat

def main():
	#Parse arguements
	parseArg(sys.argv[1:])
	
	#Open file
	f = open(inputfile,'r')
	#t = float(f.readline())
	YR2S = 3.15569e7
	#t = float(linecache.getline(inputfile, 1))
	
	#Get data
	data = np.genfromtxt(inputfile,skip_header=0)

	if data.shape[1] == 7:
		x = data[:,0]
		y = data[:,1]
		den = data[:,2]
		pre = data[:,3]
		hii = data[:,4]
		u = data[:,5]
		v = data[:,6]

		# Set up a regular grid of interpolation points
		xi, yi = np.linspace(x.min(), x.max(), 300), np.linspace(y.min(), y.max(), 300)
		xi, yi = np.meshgrid(xi, yi)

		# Interpolate; method={'nearest','linear','cubic'}
		deni = scipy.interpolate.griddata((x, y), den, (xi, yi), method='nearest')
		prei = scipy.interpolate.griddata((x, y), pre, (xi, yi), method='nearest')
		hiii = scipy.interpolate.griddata((x, y), hii, (xi, yi), method='nearest')

		ui = scipy.interpolate.griddata((x, y), u, (xi, yi), method='nearest')
		vi = scipy.interpolate.griddata((x, y), v, (xi, yi), method='nearest')

		#print np.max(ui)
		#print np.max(vi)

		#plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
		#rho
		#plt.imshow(zi, vmin=0.1, vmax=1.8, origin='lower', extent=[x.min(), x.max(), y.min(), y.max()])
		#pg
		plt.imshow(hiii,vmin=hii.min(),vmax=hii.max(),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()])
		
		#im = plt.pcolormesh(xi,yi,zi, cmap='hot', norm=LogNorm())	
		#plt.colorbar(im,orientation='vertical',format=LogFormatterMathtext())

		# Quiver Plot for vector field
		plt.quiver(xi[::9,::9], yi[::9,::9], ui[::9,::9], vi[::9,::9], pivot='mid', units='inches', color='r', scale=20)
	
	elif data.shape[1] == 5:
		x = data[:,0]
		den = data[:,1]
		pre = data[:,2]
		hii = data[:,3]
		u = data[:,4]
		#plt.xlim([0,128])
		#plt.ylim([307411,450000]) #pg
		#plt.yscale('log')
		plt.plot(x,den)
		#sod_x = [0] * 6000
		#sod_rho = [0] * 6000
		#sod_x, sod_rho = sod.sod()
		#for i in range(0, 6000):
		#	sod_x[i] = 100*sod_x[i]
		#plt.plot(sod_x, sod_rho)
	else:
		print "Can not handle %d columns" % data.shape[1]
		return
	
	#title = "%5.3f * trec" % t
	#plt.title(title)
	plt.savefig(outputfile, format=figformat)
	print 'Plotted ' + outputfile

main()


