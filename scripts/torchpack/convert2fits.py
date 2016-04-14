import numpy as np
from astropy.io import fits
import os
import argparse

import torch

###	Parse arguements
parser = argparse.ArgumentParser(description='Converts Torch raw data type to FITS image.')
parser.add_argument('inputfilename', metavar='inputfilename', type=str, help='File to convert.')
args = parser.parse_args()

raw_data = torch.CFD_Data(args.inputfilename, axial=False)

hdu = fits.PrimaryHDU()
hdu.data = raw_data.interpolate(raw_data.get_var_raw("nh"), 'nearest')

outputfilename = os.path.splitext(args.inputfilename)[0] + ".fits"
hdu.writeto(outputfilename, clobber=True)