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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import sod

inputfile = ''
outputfile = ''
paramsfile = ''
figformat = 'png'
var_type = 'pre'
quiverOn = False
DPI = 300
zlog = True

font = {'family':'serif','size':16}
plt.rc('font',**font)
plt.rc('legend',**{'fontsize':14})
plt.rc('xtick',**{'labelsize':12})
plt.rc('ytick',**{'labelsize':12})

def make_colormap(seq):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).
    """
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

def parseArg(argv, figformat):
	global inputfile
	global outputfile
	global var_type
	try:
		opts, args = getopt.getopt(argv,"i:phdt")
	except getopt.GetoptError:
		print 'ERROR: plot.py -i <inputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-i"):
			inputfile = arg
		if opt in ("-h"):
			var_type = "hii"
		if opt in ("-p"):
			var_type = "pre"
		if opt in ("-d"):
			var_type = "den"
		if opt in ("-t"):
			var_type = "tem"
	if inputfile == '':
		print 'plot.py -i <inputfile>'
		sys.exit(2)
	elif not '.' in inputfile:
		inputfile += '.xq'
	outputfile = inputfile.partition('.')[0] + '.' + figformat

def main():
	fig=plt.figure()
	#fig.set_size_inches(10,10)
	#ax=fig.add_subplot(1,1,1)

	fb0 = 0,0,0
	fb1 = 21/255.0,4/255.0,168/255.0
	fb2 = 0,247/255.0,214/255.0
	fb = make_colormap([fb0, fb1, 0.5, fb1, fb2])

	#Parse arguements
	parseArg(sys.argv[1:], figformat)
	
	#Open file
	f = open(inputfile,'r')
	#t = float(f.readline())
	YR2S = 3.15569e7
	#t = float(linecache.getline(inputfile, 1))
	
	#Get data
	data = np.genfromtxt(inputfile,skip_header=0)

	if data.shape[1] == 7:
		# Data set up.
		x = data[:,0]
		y = data[:,1]
		if var_type == 'den':
			var = data[:,2]
		if var_type == 'pre':
			var = data[:,3]
		if var_type == 'hii':
			var = data[:,4]
		if var_type == 'tem':
			var = (data[:,3]/data[:,2])*1.91939e10/(8.314462e7*(data[:,4] + 1))
		u = data[:,5]
		v = data[:,6]

		# Interpolation set up.
		xi, yi = np.linspace(x.min(), x.max(), 1024), np.linspace(y.min(), y.max(), 2048)
		xi, yi = np.meshgrid(xi, yi)
		interp = 'nearest' #{'nearest','linear','cubic'}
		vari = scipy.interpolate.griddata((x, y), var, (xi, yi), method=interp)
		ui = scipy.interpolate.griddata((x, y), u, (xi, yi), method=interp)
		vi = scipy.interpolate.griddata((x, y), v, (xi, yi), method=interp)

		# Figure set up.
		ax = fig.add_axes([0.04, 0.04, 0.92, 0.92])
		divider = make_axes_locatable(ax)
		c_ax = divider.append_axes("right", size="5%", pad=0.04)
		vari_min = vari.min()
		vari_max = vari.max()
		im = ax.imshow(vari,vmin=var.min(),vmax=var.max(),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()], interpolation="none", cmap=fb)
		cbar = plt.colorbar(im, cax=c_ax)

		# Quiver plot for vector field
		if quiverOn == True:
			plt.quiver(xi[::9,::9], yi[::9,::9], ui[::9,::9], vi[::9,::9], pivot='mid', units='inches', color='r', scale=0.2)
	
	elif data.shape[1] == 5:
		x = data[:,0]
		if var_type == 'den':
			var = data[:,1]
		if var_type == 'pre':
			var = data[:,2]
		if var_type == 'hii':
			var = data[:,3]
		if var_type == 'tem':
			var = (data[:,2]/data[:,1])*1.91939e10/(8.314462e7*(data[:,3] + 1))
		u = data[:,4]
		#plt.xlim([0,128])
		#plt.ylim([307411,450000]) #pg
		#plt.yscale('log')
		plt.ylim([0, max(var)]);
		plt.plot(x,var)
		#sod_x = [0] * 6000
		#sod_rho = [0] * 6000
		#sod_x, sod_rho = sod.sod()
		#for i in range(0, 6000):
		#	sod_x[i] = 100*sod_x[i]
		#plt.plot(sod_x, sod_rho)
	else:
		print "Can not handle %d columns" % data.shape[1]
		return

	#extent = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
	#fig.savefig(outputfile, format=figformat, dpi=DPI, bbox_inches=extent, pad_inches=0)
	fig.savefig(outputfile, format=figformat, dpi=DPI, bbox_inches='tight', pad_inches=0.1)
	print 'plot.py: plotted ' + var_type + ' in ' + outputfile

main()


