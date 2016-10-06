import numpy as np
import lupa
lua = lupa.LuaRuntime(unpack_returned_tuples=True)

# Parse arguments.
import argparse
parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
args = parser.parse_args()

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("cdat", "../torchpack/cornishdata.py")
import cdat
cornish_data = cdat.CornishData(args.iden)
star_data = cornish_data.star_data
nstars = len(star_data[:,0])

data = np.genfromtxt(cornish_data.dirname + "/total-fluxes.txt", skip_header=1)

remaining_set = np.zeros(nstars + 1, dtype=bool)

sigma = 7

for i in range(0, len(data[:,0])):
	index = int(data[i,0])
	tot_flux = data[i,1]

	ipad = cornish_data.star_id_str(index)
	dirname = cornish_data.dirname + "/star_" + ipad + "/"

	filestring = open(dirname + "radioconfig_0.lua", 'r').read()
	table = lua.eval("{" + filestring + "}")
	dec = table["Parameters"]["declination"]

	if (dec >= 14.2 and tot_flux >= 0.25 * sigma) or (dec < 14.2 and tot_flux >= 0.35 * sigma):
		remaining_set[index] = True

ofile1 = open(cornish_data.dirname + "/excl-weak-fluxes.txt", 'w')
ofile1.write("star_id [(dec >= 14.2 deg && tot_flux < 5 * 0.25) || (dec < 14.2 deg && tot_flux < 5 * 0.35)]\n")

for index in range(1, nstars + 1):
	if not remaining_set[index]:
		ofile1.write(str(int(index)) + "\n")

ofile1.close()
