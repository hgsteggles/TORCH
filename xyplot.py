#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import sys, getopt
import math
import linecache
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext

paramsfile = ''
figformat = 'png'
zlog = True

x = data[:,0]
y = data[:,1]
den = data[:,2]

# Set up a regular grid of interpolation points
xi, yi = np.linspace(x.min(), x.max(), 300), np.linspace(y.min(), y.max(), 300)
xi, yi = np.meshgrid(xi, yi)

# Interpolate; method={'nearest','linear','cubic'}
deni = scipy.interpolate.griddata((x, y), den, (xi, yi), method='nearest')
plt.imshow(vari,vmin=var.min(),vmax=var.max(),origin='lower',extent=[x.min(),x.max(),y.min(),y.max()])

plt.savefig(outputfile, format=figformat, interpolate=none)
print 'Plotted ' + outputfile


