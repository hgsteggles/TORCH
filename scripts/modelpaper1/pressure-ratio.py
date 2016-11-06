import numpy as np
from astropy.io import fits
import lupa

lua = lupa.LuaRuntime(unpack_returned_tuples=True)

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
import torch

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData(useSet4=True)

itime = 25

def ratio_mvd():
	for index in range(1, 46):
		col = int(index / 9)
		row = index % 9

		filename = mp1_data.getDataFilenameByIndex(index, itime)
		datacube = torch.CFD_Data(filename, axial=False)

		vol = datacube.get_var_raw('vol-cylindrical')
		pressure = datacube.get_var_raw('pre') * vol
		tem = datacube.get_var_raw('tem')

		wind_indexes = np.where(tem > 12000)
		hii_indexes = np.where((tem > 7000) & (tem < 12000))

		wind_vol = np.sum(vol[wind_indexes])
		hii_vol = np.sum(vol[hii_indexes])

		wind_pressure = np.sum(pressure[wind_indexes]) / wind_vol
		hii_pressure = np.sum(pressure[hii_indexes]) / hii_vol

		ratio = wind_pressure / hii_pressure

		print str(index) + " " +  str(ratio)


def ratio_mvt():
	snapshots = [10, 20, 30, 40, 50]

	for col in range(5):
		for row in range(9):
			filename = mp1_data.getDataFilename(row, 2, snapshots[col])
			datacube = torch.CFD_Data(filename, axial=False)

			vol = datacube.get_var_raw('vol-cylindrical')
			pressure = datacube.get_var_raw('pre') * vol
			tem = datacube.get_var_raw('tem')

			wind_indexes = np.where(tem > 12000)
			hii_indexes = np.where((tem > 7000) & (tem < 12000))

			wind_vol = np.sum(vol[wind_indexes])
			hii_vol = np.sum(vol[hii_indexes])

			wind_pressure = np.sum(pressure[wind_indexes]) / wind_vol
			hii_pressure = np.sum(pressure[hii_indexes]) / hii_vol

			ratio = wind_pressure / hii_pressure

			print str(col) + " " + str(row) + " "  +  str(ratio)

ratio_mvd()