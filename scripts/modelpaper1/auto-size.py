import numpy as np
from astropy.io import fits
import math

import lupa
lua = lupa.LuaRuntime(unpack_returned_tuples=True)

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData()

fileprefix2 = "data/model_paper1/set2/"
fileprefix3 = "data/model_paper1/set3/"

offgrid = np.fromfile(fileprefix2 + "offgrid_times", dtype=int, sep='\n')

def getFilename(i, set2):
	padi = "%02d" % ((i+1),)
	filesuffix = "data_" + padi + "/data2D_025.fits"

	if set2:
		return fileprefix2 + filesuffix
	else:
		return fileprefix3 + filesuffix

def getParamFilename(i, set2):
	padi = "%02d" % ((i+1),)
	filesuffix = "data_" + padi + "/parameters_" + str(i+1) + ".lua"

	if set2:
		return fileprefix2 + filesuffix
	else:
		return fileprefix3 + filesuffix

def calc_sizes(index, itime):
	filename = mp1_data.getDataFilenameByIndex(index, itime)
	paramname = mp1_data.getParamFilenameByIndex(index,itime)

	filestring = open(paramname, 'r').read()
	table = lua.eval("{" + filestring + "}")
	L = table["Parameters"]["Grid"]["side_length"]
	nx = table["Parameters"]["Grid"]["no_cells_x"]
	dx = L / float(nx)

	datacube = torch.CFD_Data(filename, axial=False)
	tem = datacube.get_var_raw('tem')

	nhii = (tem > 4000).sum() * 2
	area_hii = nhii * dx * dx
	radius_hii = np.sqrt(area_hii / math.pi)

	nwind = (tem > 20000).sum() * 2
	area_wind = nwind * dx * dx
	radius_wind = np.sqrt(area_wind / math.pi)

	return radius_wind, radius_hii

rad_sw, rad_if = calc_sizes(24, 30)

print rad_if / 3.09e18
print rad_sw / 3.09e18