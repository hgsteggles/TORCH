import argparse
import os

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")

import torch


###	Parse arguements
parser = argparse.ArgumentParser(description='Converts Torch data to FITS format.')
parser.add_argument('-s', action='store_true', help='Run script silently.')
parser.add_argument('var_type', metavar='var_type', type=str, choices=torch.CFD_Data.var_typenames, help='Variable type.')
parser.add_argument('inputfile', metavar='inputfile', type=str, help='Filename of input Torch data.')
parser.add_argument('-o', metavar='outputfile', type=str, help='Filename of output FITS data.')
args = parser.parse_args()

data = torch.CFD_Data(args.inputfile)

outputfile = args.o
if outputfile == None:
    outputfile = os.path.splitext(args.inputfile)[0] + ".fits"
    print outputfile

data.convert_to_fits(torch.VarType(args.var_type, False), outputfile)