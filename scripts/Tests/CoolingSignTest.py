import numpy as np
import sys
import glob
import argparse
import math
import linecache
from joblib import Parallel, delayed

from torch_plot import TorchCFD, TorchCool

parser = argparse.ArgumentParser(description='Tests validity of cooling data.')

parser.add_argument('coolfile', metavar='coolfile', type=str, help='Cooling file.')
args = parser.parse_args()

cool = TorchCool(args.coolfile)

def test_process(vname):
	global cool
	v = cool.get_var(vname)
	if len(v[v > 0]) > 0:
		print vname + " has positive values!"

test_process("imlc")
test_process("nmlc")
test_process("recc")
test_process("colc")
test_process("ciec")
test_process("nmc")
