import numpy as np
from astropy.io import fits
import lupa

lua = lupa.LuaRuntime(unpack_returned_tuples=True)

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

ofile = open("central-pressures", 'w')

for index in range(45):
	isSet2 = offgrid[index] > 25
	filename = getFilename(index, isSet2)

	filestring = open(getParamFilename(index, isSet2), 'r').read()
	table = lua.eval("{" + filestring + "}")
	xs = table["Parameters"]["Star"]["cell_position_x"]
	ys = table["Parameters"]["Star"]["cell_position_y"]

	hdu_list = fits.open(filename)
	image_data = hdu_list[0].data

	nh = image_data[ys][xs]

	ofile.write(str(index + 1) + '\t' + str(nh) + '\n')