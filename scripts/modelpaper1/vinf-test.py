import numpy as np
import math
import sys

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("koo", "../torchpack/koo.py")
load_src("henney", "../misc/henney-cooling.py")

import koo
import henney

load_src("torch", "../torchpack/torch.py")
import torch

import lupa
lua = lupa.LuaRuntime(unpack_returned_tuples=True)

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData(useSet4=True)

CM2PC = 3.24078e-19
PC2CM = 1.0/CM2PC
R_GAS = 8.3144598e7
mh = 1.6738232e-24
S2YR = 3.17098e-8
YR2S = 3.154e7
KB = 1.3806485e-16

fileprefix2 = "data/model_paper1/set2/"
fileprefix3 = "data/model_paper1/set3/"

snapshots = [10, 20, 30, 40, 50]

vinf = []

for i in range(45):
	col = int(i / 9)
	row = i % 9
	itime = 50 if i >= 19 and i <= 27 else 20

	filename = mp1_data.getParamFilename(row, col, itime)
	filestring = open(filename, 'r').read()

	table = lua.eval("{" + filestring + "}")
	vinf.append(table["Parameters"]["Star"]["wind_velocity"])

itime = 50
for index in range(1, 45):
	filename = mp1_data.getDataFilenameByIndex(index, itime)

	datacube = torch.CFD_Data(mp1_data.getDataFilenameByIndex(index, itime), axial=False)
	maxvtot = np.max(datacube.get_var_raw('vtot'))

	rel = abs(100.0 * (maxvtot - float(vinf[index - 1])) / float(vinf[index - 1]))

	print str(index) + " " + str(maxvtot / vinf[index - 1])

