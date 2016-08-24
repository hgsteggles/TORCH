import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
import torch
load_src("hgspy", "../torchpack/hgspy.py")
import hgspy
load_src("mp1", "../torchpack/modelpaper1.py")
import mp1

modeldata = mp1.ModelData()

masses = np.array(modeldata.masses, dtype=float)
mvt = np.array(modeldata.mvt_if_r, dtype=float)
y = []

for i in range(len(modeldata.mvt_if_r)):
	slope, intercept, r_value, p_value, std_err = stats.linregress(masses, mvt[i])

	def func(x):
		return slope * x + intercept

	print std_err

	y.append(func(masses))

plt.plot(modeldata.masses, y[0])

plt.show()
