import linecache
import numpy as np
import scipy
import matplotlib.pyplot as plt

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("torch", "../torchpack/torch.py")
load_src("hgspy", "../torchpack/hgspy.py")

import torch
import hgspy

def meshgrid2(*arrs):
    arrs = tuple(reversed(arrs))  #edit
    lens = map(len, arrs)
    dim = len(arrs)

    sz = 1
    for s in lens:
        sz*=s

    ans = []
    for i, arr in enumerate(arrs):
        slc = [1]*dim
        slc[i] = lens[i]
        arr2 = np.asarray(arr).reshape(slc)
        for j, sz in enumerate(lens):
            if j != i:
                arr2 = arr2.repeat(sz, axis=j)
        ans.append(arr2)

    return tuple(ans[::-1])

DPI = 300
figformat = 'png'
plot_size = 5
fontsize = 16
torch.set_font_sizes(fontsize=16)

inputfile = "data/galsim/density3D.txt"

### Read data.
dimstr = linecache.getline(inputfile, 1).split()
nx = int(dimstr[0])
ny = int(dimstr[1])
nz = int(dimstr[2])

data = np.genfromtxt(inputfile, skip_header=2)

den = data[:,5]

x = []
x.append(data[:,0])
x.append(data[:,1])
x.append(data[:,2])

dx = []
dx.append(x[0][1] - x[0][0])
dx.append(x[1][1] - x[1][0])
dx.append(x[2][1] - x[2][0])

xi = []
xi.append(np.linspace(x[0].min(), x[0].max(), nx))
xi.append(np.linspace(x[1].min(), x[1].max(), ny))
xi.append(np.linspace(x[2].min(), x[2].max(), nz))
xi[0], xi[1], xi[2] = meshgrid2(xi[0], xi[1], xi[2])

deni = scipy.interpolate.griddata((x[0], x[1], x[2]), den,
						   (xi[0], xi[1], xi[2]),
						   method="nearest")
coli = np.zeros(shape=deni[:,:,0].shape)
for i in range(len(deni[0,0,:])):
	coli += deni[:,:,i]

def plot_image(var):
	return plt.imshow(var, vmin=var.min(), vmax=var.max(), origin='lower',
			   extent=[x[0].min(),x[0].max(),x[1].min(),x[1].max()],
			   interpolation="nearest", cmap=hgspy.get_par_cmap())

im = plot_image(deni[:,:,15])

plt.show()

### Plotting
#plotter = torch.Plotter(nx, ny, plot_size, figformat, DPI)
#plotparams = torch.PlotParams(None, None, None, False, 'linear', (1, 1), None, tight=False, detail="all")
#grid = plotter.getGrid(plotparams)
#plotter.modifyGrid(grid, True)

#ax = grid[0]
