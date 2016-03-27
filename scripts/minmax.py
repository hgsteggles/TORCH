import numpy as np
import sys
import os
import glob
import argparse
import math
import linecache
from joblib import Parallel, delayed

from torch_plot import TorchCFD, TorchCool

def safe_log10(v):
	if v <= 0:
		return -float("inf")
	elif not np.isfinite(v):
		return v
	else:
		return math.log10(v)

parser = argparse.ArgumentParser(description='Finds variable range of all CFD cooling data in directory.')

parser.add_argument('-cfd', metavar='cfdfile', type=str, help='List of cfd files.')
args = parser.parse_args()

vmina = []
vmaxa = []

filenames = glob.glob(args.cfd)

def minmax(cfdfile):
	###	Get data
	cfd = TorchCFD(cfdfile)

	nh = cfd.get_var_raw('nh')
	pre = cfd.get_var_raw('pre')
	nhii = cfd.get_var_raw('nhii')
	tem = cfd.get_var_raw('tem')

	return [nh.min(), pre.min(), nhii.min(), tem.min(), nh.max(), pre.max(), nhii.max(), tem.max()]

minmaxes = np.array(Parallel(n_jobs=8)(delayed(minmax)(f) for f in filenames))

mins = np.amin(minmaxes[:,:4], axis=0)
maxs = np.amax(minmaxes[:,4:], axis=0)

ofilename = os.path.dirname(os.path.realpath(filenames[0])) + '/minmax.log'
ofile = open(ofilename, 'w')

ofile.write("nh         = [ " + str(mins[0]) + "\t" + str(maxs[0]) + " ]\n")
ofile.write("pre        = [ " + str(mins[1]) + "\t" + str(maxs[1]) + " ]\n")
ofile.write("nhii       = [ " + str(mins[2]) + "\t" + str(maxs[2]) + " ]\n")
ofile.write("tem        = [ " + str(mins[3]) + "\t" + str(maxs[3]) + " ]\n\n")

ofile.write("log10_nh   = [ " + str(safe_log10(mins[0])) + "\t" + str(safe_log10(maxs[0])) + " ]\n")
ofile.write("log10_pre  = [ " + str(safe_log10(mins[1])) + "\t" + str(safe_log10(maxs[1])) + " ]\n")
ofile.write("log10_nhii = [ " + str(safe_log10(mins[2])) + "\t" + str(safe_log10(maxs[2])) + " ]\n")
ofile.write("log10_tem  = [ " + str(safe_log10(mins[3])) + "\t" + str(safe_log10(maxs[3])) + " ]\n")


