import argparse
import os
import sys
import matplotlib.pyplot as plt
import numpy as np

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
import torch

###	Parse arguements
parser = argparse.ArgumentParser(description='Plots 2D image of CFD data.')
parser.add_argument('inputfile', metavar='inputfile', type=str, help='Input file to find integrate variable over.')
parser.add_argument('var_type', metavar='var_type', type=str, choices=torch.Cool_Data.var_typenames, help='Variable type.')
args = parser.parse_args()

###	Data set up
datacube = torch.Cool_Data(args.inputfile, axial=False)

dx =  ((0.5 * 3.09e18) / 150)
vol = datacube.get_var_raw('vol-cylindrical') * (3.09e18)**3 * dx**3
voltot =  np.sum(vol)

print np.sum(datacube.get_var_raw(args.var_type) * vol)


