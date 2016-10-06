import linecache
import numpy as np
import scipy
import scipy.stats
import matplotlib.pyplot as plt
import math

# Parse arguments.
import argparse

parser = argparse.ArgumentParser(description='Calculates max flux in cornish set.')
parser.add_argument('iden', metavar='iden', type=int, help='Density index in CORNISH param set.')
args = parser.parse_args()


def load_src(name, fpath):
	import os, imp
	return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))


load_src("torch", "../torchpack/torch.py")
import torch

load_src("hgspy", "../torchpack/hgspy.py")
import hgspy

load_src("koo", "../torchpack/koo.py")
import koo

DPI = 300
figformat = 'png'
plot_size = 5
fontsize = 16
torch.set_font_sizes(fontsize=16)

YR2S = 3.154e7
PC2CM = 3.086e18
CM2PC = 1.0 / PC2CM

inputfile = "data/galsim/starpopfinal.txt"

high_mass = True
windowed = False
age_cutoff = True
size_cut = False
flux_cut = False

suffix = ""
if high_mass:
	suffix += "-hm"
if windowed:
	suffix += "-w"
if age_cutoff:
	suffix += "-a"
if size_cut:
	suffix += "-s"
if flux_cut:
	suffix += "-f"
outputfile = "data/galsim/starpop" + suffix + ".txt"

firstline = linecache.getline(inputfile, 1)

###	Data set up.
params = np.genfromtxt("config/zams-stripped.txt", skip_header=1)

massdat = params[:, 0]  # [Msun]
logQdat = params[:, 6]  # log[Q / phot.s-1]

data = np.genfromtxt(inputfile, skip_header=1)

age_index = 2

# Select high mass stars (min mass = 6 Msun).
if high_mass:
	mcur = data[:, 0]
	data = np.delete(data, np.where(mcur < 6), 0)
# Select stars in survey window 10 < l < 65 and |b| < 1
if windowed:
	gall = data[:, 15]
	data = np.delete(data, np.where(gall < 10), 0)
	gall = data[:, 15]
	data = np.delete(data, np.where(gall > 65), 0)
	galb = data[:, 16]
	data = np.delete(data, np.where(np.abs(galb) > 1), 0)
if age_cutoff:
	age = data[:, age_index]
	data = np.delete(data, np.where(age > 800), 0)
	age = data[:, age_index]
	data = np.delete(data, np.where(age <= 0), 0)
if size_cut:
	mcur = data[:, 0]
	dsun = data[:, 4] * 1000.0 * PC2CM
	age = data[:, age_index] * 1000.0 * YR2S
	nh = 0.8e4
	logQ = scipy.interpolate.interp1d(massdat, logQdat)(mcur)
	raga = koo.calcSpitzerRadius(logQ, nh, age)
	stag = koo.calcStagnationRadius(logQ, nh)
	radius = np.minimum(raga, stag)
	angsize = (radius / dsun) * (180.0 / math.pi) * 60.0 * 60.0

	data = np.delete(data, np.where(0.75 * angsize > 15.0), 0)


def calcPlanck(tem, freq):
	h = 6.626e-34
	c = 3.0e8
	kb = 1.38e-23
	return (2.0 * h * freq * freq * freq / (c * c)) / (math.exp(h * freq / (kb * tem)) - 1)


def calcTau(T, freq, ne, logq):
	return 0.53 * (T / 1.0e4) ** (-1.35) * (freq / 1.0e9) ** (-2.1) * (ne / 1.0e3) ** 0.66 * np.power(10.0, (logq - 49.7) / 3.0)


nhs = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]

if flux_cut:
	mcur = data[:, 0]
	dsun = data[:, 4] * 1000.0 * PC2CM
	age = data[:, age_index] * 1000.0 * YR2S
	nh = nhs[args.iden]
	logQ = scipy.interpolate.interp1d(massdat, logQdat)(mcur)
	stro = koo.calcSpitzerRadius(logQ, nh, age)
	stag = koo.calcStagnationRadius(logQ, nh)
	radius = np.minimum(stro, stag)
	angsize_rad = (radius / dsun)
	angsize_as = angsize_rad * 60.0 * 60.0 * 180.0 / math.pi
	ne = nh * (stro / radius) ** (1.5)
	T = 8000.0
	f = 5.0e9
	tau = calcTau(T, f, ne, logQ)
	B = calcPlanck(T, f)
	S = B * (1.0 - np.exp(-tau))
	F = S * math.pi * angsize_rad * angsize_rad * 1.0e26 * 1000

	print F

	data = np.delete(data, np.where(20 * F < 5 * 0.25), 0)

file = open(outputfile, 'w')
file.write(firstline)
np.savetxt(file, data, fmt='%-9.9s')
