import numpy as np
import sys
import glob
import argparse
import math
import linecache
from joblib import Parallel, delayed

from torch_plot import TorchCFD, TorchCool

parser = argparse.ArgumentParser(description='Finds variable range of all CFD cooling data in directory.')

parser.add_argument('-cfd', metavar='cfdfile', type=str, help='List of cfd files.')
parser.add_argument('-cool', metavar='coolfile', type=str, help='List of cooling files.')
args = parser.parse_args()

vmina = []
vmaxa = []

filenames = np.array([glob.glob(args.cfd), glob.glob(args.cool)])
filenames = filenames.transpose()

def minmax(cfdfile, coolfile):
	###	Get data
	cfd = TorchCFD(cfdfile)
	cool = TorchCool(coolfile)

	nh = cfd.get_var('nh')
	tem = cfd.get_var('tem')
	crate = cool.get_var('cool')
	L = np.absolute(cool.get_var('lcool')/(nh*nh))
	E = 1.5*nh*1.3806488e-16*tem
	with np.errstate(divide='ignore', over='ignore'):
		dt = E / np.absolute(crate) / 3.15569e7
		dt[crate == 0] = 1e10

	return [nh.min(), tem.min(), E.min(), L.min(), dt.min(), nh.max(), tem.max(), E.min(), L.max(), dt.max()]

minmaxes = np.array(Parallel(n_jobs=8)(delayed(minmax)(f[0], f[1]) for f in filenames))

mins = np.amin(minmaxes[:,:5], axis=0)
maxs = np.amax(minmaxes[:,5:], axis=0)

print "nh        = [ " + str(mins[0]) + "\t" + str(maxs[0]) + " ]"
print "tem       = [ " + str(mins[1]) + "\t" + str(maxs[1]) + " ]"
print "E         = [ " + str(mins[2]) + "\t" + str(maxs[2]) + " ]"
print "L         = [ " + str(mins[3]) + "\t" + str(maxs[3]) + " ]"
print "dt        = [ " + str(mins[4]) + "\t" + str(maxs[4]) + " ]"

def safe_log10(v):
	if v <= 0:
		return -float("inf")
	elif not np.isfinite(v):
		return v
	else:
		return math.log10(v)

print "log10_nh  = [ " + str(safe_log10(mins[0])) + "\t" + str(safe_log10(maxs[0])) + " ]"
print "log10_tem = [ " + str(safe_log10(mins[1])) + "\t" + str(safe_log10(maxs[1])) + " ]"
print "log10_E   = [ " + str(safe_log10(mins[2])) + "\t" + str(safe_log10(maxs[2])) + " ]"
print "log10_L   = [ " + str(safe_log10(mins[3])) + "\t" + str(safe_log10(maxs[3])) + " ]"
print "log10_dt  = [ " + str(safe_log10(mins[4])) + "\t" + str(safe_log10(maxs[4])) + " ]"

#print "dt_cool:     " + str(vmin) + " " + str(vmax)
#print "log_dt_cool: " + str(math.log10(vmin)) + " " + str(math.log10(vmax))
