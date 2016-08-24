import math
from scipy import optimize

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
import torch
load_src("hgspy", "../torchpack/hgspy.py")
import hgspy
load_src("koo", "../torchpack/koo.py")
import koo

CM2PC = 3.24078e-19
PC2CM = 1.0 / CM2PC
YR2S = 3.154e7

times = [2.0e4, 4.0e4, 6.0e4, 8.0e4, 10.0e4]


def find_stromgren(logQ, nh):
	alpha = 2
	rstar = 0.35 * PC2CM
	rc = 0.01 * PC2CM
	a = 1 + (rstar / rc)**2.0
	n0 = nh * a**(alpha / 2.0)

	a = 1 + (rstar / rc)**2.0
	b = (10.0**logQ) / (4.0 * math.pi * koo.calcAlphaB(8000.0) * (rc**3) * (n0**2))

	def func(x):
		return -x / (2.0 * (a + x**2)) + math.atan2(x, math.sqrt(a)) / (2.0 * math.sqrt(a)) - b

	return optimize.brentq(func, 0, 100) * rc * CM2PC

RS = find_stromgren(49.81, 3.2e4)
print RS
print koo.calcSpitzerRadius2(RS * PC2CM, times[0]*YR2S) * CM2PC