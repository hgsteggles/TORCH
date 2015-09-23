import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import linecache

from mpl_toolkits.axes_grid1 import AxesGrid

class TorchData:
	var_types = []

	def __init__(self, inputfile, axial=False):
		self.t = int(float(linecache.getline(inputfile, 1))/31536000.0/10)*10
		self.nx = (2 if axial else 1)*int(linecache.getline(inputfile, 2))
		self.ny = int(linecache.getline(inputfile, 3))
		self.nz = int(linecache.getline(inputfile, 4))
		self.data = np.genfromtxt(inputfile, skip_header=4)

		self.nd = 3
		if self.nz == 1:
			self.nd = 2
		if self.ny == 1:
			self.nd = 1

		if (self.nd == 3):
			minZ = z.min()
			self.data = data[z == minZ]

		if axial:
			self.mirror_data()

		self.x = []
		coords = ['x', 'y']
		for i in range(min(2, self.nd)):
			self.x.append(self.get_var(coords[i]))

		self.xi = []
		self.xi.append(np.linspace(self.x[0].min(), self.x[0].max(), self.nx))
		if self.nd > 1:
			self.xi.append(np.linspace(self.x[1].min(), self.x[1].max(), self.ny))
			self.xi[0], self.xi[1] = np.meshgrid(self.xi[0], self.xi[1])

	def mirror_data(self):
		data2 = self.data[self.get_var('x') != 0]
		data2[:,0] = -data2[:,0]
		self.data = np.append(self.data, data2, axis=0)

	def get_var(self, var_type):
		print "TorchData::get_var: " + var_type + " does not exist in this data."
		return None

	def safe_log10(self, var):
		try:
			var[var == 0] = var[var != 0].min()
		except ValueError:
			pass
		return np.log10(var)

	def get_log10_variable(self, var_name):
		return safe_log10(get_variable(self, var_name))

	def appropriate_to_log(self, var_type):
		return True

	def interpolate(self, var, interp):
		return scipy.interpolate.griddata((self.x[0], self.x[1]), var, (self.xi[0], self.xi[1]), method=interp)

class TorchCFD(TorchData):
	var_types = ['x', 'y', 'z', 'nh', 'den', 'pre', 'hii', 'nhii', 'nhi', 'tem', 'vel0', 'vel1']

	def __init__(self, inputfile, axial=False):
		TorchData.__init__(self, inputfile, axial)

	def mirror_data(self):
		data2 = self.data[self.get_var('x') != 0]
		data2[:,0] = -data2[:,0]
		data2[:,self.nd+3] = -data2[:,self.nd+3]
		self.data = np.append(self.data, data2, axis=0)
		
	def get_var(self, var_type):
		if self.nd != 0:
			if var_type == 'x':
				return self.data[:,0]/3.09e18
			elif var_type == 'y' and self.nd > 1:
				return self.data[:,1]/3.09e18
			elif var_type == 'z' and self.nd > 2:
				return self.data[:,2]/3.09e18
			elif var_type == 'nh':
				return self.get_var('den')/1.674e-24
			elif var_type == 'den':
				return self.data[:,self.nd]
			elif var_type == 'pre':
				return self.data[:,self.nd+1]
			elif var_type == 'hii':
				return self.data[:,self.nd+2]
			elif var_type == 'nhii':
				return self.get_var('hii')*self.get_var('den')/1.674e-24
				#return self.data[:,self.nd+2]*self.data[:,self.nd]/1.674e-24
			elif var_type == 'tem':
				return (self.get_var('pre')/self.get_var('den'))/(8.314462e7*(self.get_var('hii') + 1))
				#return (self.data[:,self.nd+1]/self.data[:,self.nd])/(8.314462e7*(self.data[:,self.nd+2] + 1))
			elif var_type == 'vel0':
				return self.data[:,self.nd+3]
			elif var_type == 'vel1' and self.nd > 1:
				return self.data[:,self.nd+4]
			elif var_type == 'vel2' and self.nd > 2:
				return self.data[:,self.nd+5]
			else:
				return TorchData.get_var(self, var_type)

	def appropriate_to_log(self, var_type):
		return var_type != 'hii' and var_type != 'vel0' and var_type != 'vel1' and var_type != 'vel2'

class TorchCool(TorchData):
	var_types = ['x', 'y', 'z', 'imlc', 'nmlc', 'recc', 'colc', 'ciec', 'nmc', 'fuvh', 'irh', 'crh', 'heat', 'cool', 'tot']

	def __init__(self, inputfile, axial=False):
		TorchData.__init__(self, inputfile, axial)

	def mirror_data(self):
		data2 = self.data[self.get_var('x') != 0]
		data2[:,0] = -data2[:,0]
		self.data = np.append(self.data, data2, axis=0)
		
	def get_var(self, var_type):
		if self.nd != 0:
			if var_type == 'x':
				return self.data[:,0]/3.09e18
			elif var_type == 'y' and self.nd > 1:
				return self.data[:,1]/3.09e18
			elif var_type == 'z' and self.nd > 2:
				return self.data[:,2]/3.09e18
			elif var_type == 'imlc':
				return self.data[:,self.nd]
			elif var_type == 'nmlc':
				return self.data[:,self.nd+1]
			elif var_type == 'recc':
				return -self.data[:,self.nd+2]
			elif var_type == 'colc':
				return self.data[:,self.nd+3]
			elif var_type == 'ciec':
				return self.data[:,self.nd+4]
			elif var_type == 'nmc':
				return self.data[:,self.nd+5]
			elif var_type == 'fuvh':
				return self.data[:,self.nd+6]
			elif var_type == 'irh':
				return self.data[:,self.nd+7]
			elif var_type == 'crh':
				return self.data[:,self.nd+8]
			elif var_type == 'heat':
				return self.get_var('fuvh')+self.get_var('irh')+self.get_var('crh')
				#return self.data[:,self.nd+6]+self.data[:,self.nd+7]+self.data[:,self.nd+8]
			elif var_type == 'cool':
				return self.get_var('imlc')+self.get_var('nmlc')+self.get_var('recc')+self.get_var('colc')+self.get_var('ciec')+self.get_var('nmc')
				#return self.data[:,self.nd]+self.data[:,self.nd+1]+self.data[:,self.nd+2]+self.data[:,self.nd+3]+self.data[:,self.nd+4]+self.data[:,self.nd+5]
			elif var_type == 'lcool':
				return self.get_var('imlc')+self.get_var('recc')+self.get_var('colc')+self.get_var('ciec')+self.get_var('nmc')
				#return self.data[:,self.nd]+self.data[:,self.nd+2]+self.data[:,self.nd+3]+self.data[:,self.nd+4]+self.data[:,self.nd+5]
			elif var_type == 'tot':
				return self.data[:,self.nd+6]+self.data[:,self.nd+7]+self.data[:,self.nd+8]+self.data[:,self.nd]+self.data[:,self.nd+1]+self.data[:,self.nd+2]+self.data[:,self.nd+3]+self.data[:,self.nd+4]+self.data[:,self.nd+5]

class TorchInterpolatedGrid:
	def __init__(self, x, y, nx, ny):
		self.nx = nx
		self.ny = ny
		self.x = x
		self.y = y
		self.xi, self.yi = np.linspace(x.min(), x.max(), nx), np.linspace(y.min(), y.max(), ny)
		self.xi, self.yi = np.meshgrid(self.xi, self.yi)

	def interpolate(self, var, interp):
		return scipy.interpolate.griddata((self.x, self.y), var, (self.xi, self.yi), method=interp)

class TorchPlotter:
	def __init__(self, torch_data, plot_size=5, figformat='jpg', dpi=300):
		self.plot_size = plot_size
		self.figformat = figformat
		self.dpi = dpi
		self.cells_per_inch = torch_data.nx/plot_size
		self.aspect_ratio = torch_data.ny/float(torch_data.nx)
		self.fig = plt.figure()
		self.fig.set_size_inches(plot_size,plot_size*self.aspect_ratio)
		self.tdat = torch_data

	def image(self, var, var_min, var_max, color_map, ax):
		im = ax.imshow(var, vmin=var_min, vmax=var_max, origin='lower', extent=[self.tdat.x[0].min(),self.tdat.x[0].max(),self.tdat.x[1].min(),self.tdat.x[1].max()], interpolation="none", cmap=color_map)
		return im

	def add_quiver(self, every, ax):
		e = every
		ui = self.tdat.interpolate(self.tdat.get_var('vel0'), 'nearest')
		vi = self.tdat.interpolate(self.tdat.get_var('vel1'), 'nearest')
		max_vel = ui.max()
		max_vel = max(max_vel, vi.max())
		ax.quiver(self.tdat.xi[0][::e,::e], self.tdat.xi[1][::e,::e], ui[::e,::e], vi[::e,::e], pivot='mid', units='inches', color='r', scale=0.8*max_vel*self.cells_per_inch/e)
	
	def save_plot(self, outputfile, ax=None):
		self.fig.tight_layout()
		if ax == None:
			self.fig.savefig(outputfile, format=self.figformat, dpi=self.dpi, bbox_inches='tight', pad_inches=0.1)
		else:
			extent = ax.get_window_extent().transformed(self.fig.dpi_scale_trans.inverted())
			self.fig.savefig(outputfile, format=self.figformat, dpi=self.dpi, bbox_inches=extent, pad_inches=0)

	def single(self, vtype, var, var_min, var_max, color_map, detail="all"):
		axpad = 0.35
		cbloc = "top"
		cbmode = "each"
		if detail == "some" or detail == "none":
			axpad = 0
			cbloc = "right"
			cbmode = "single"

		grid = AxesGrid(self.fig, 111, nrows_ncols = (1, 1),
				axes_pad = axpad,
				label_mode = "L", 
				share_all = True, 
				cbar_location=cbloc,
				cbar_mode=cbmode, 
				cbar_size="6%", 
				cbar_pad="0%")
		if detail != "none":
			grid[0].set_xlabel('r / pc')
			grid[0].set_ylabel('z / pc')
		else:
			grid[0].get_xaxis().set_visible(False)
			grid[0].get_yaxis().set_visible(False)
		im = self.image(var, var_min, var_max, color_map, grid[0])
		if detail != "some" and detail != "none":
			grid.cbar_axes[0].colorbar(im)
		if detail == "some" or detail == "none":
			grid.cbar_axes[0].colorbar(im)
			grid.cbar_axes[0].toggle_label(False)
			ts = 0.001
			grid[0].text(ts, 1-ts, vtype, fontsize=10, color='red', horizontalalignment='left', verticalalignment='top', transform = grid[0].transAxes)
		if detail == "none":
			grid.cbar_axes[0].get_yaxis().set_visible(False)
		return grid[0]

	def twobytwo(self, vtypes, vs, vsminmax, color_map, detail="all"):
		axpad = 0.35
		cbloc = "top"
		cbmode = "each"
		if detail == "some" or detail == "none":
			axpad = 0
			cbloc = "right"
			cbmode = "single"

		grid = AxesGrid(self.fig, 111, nrows_ncols = (2, 2),
				axes_pad = axpad,
				label_mode = "L", 
				share_all = True, 
				cbar_location=cbloc,
				cbar_mode=cbmode, 
				cbar_size="6%", 
				cbar_pad="0%")
		for i in range(4):
			if detail != "none":
				grid[i].set_xlabel('r / pc')
				grid[i].set_ylabel('z / pc')
			else:
				grid[i].get_xaxis().set_visible(False)
				grid[i].get_yaxis().set_visible(False)
			im = self.image(vs[i], vsminmax[i][0], vsminmax[i][1], color_map, grid[i])
			if detail != "some" and detail != "none":
				grid.cbar_axes[i].colorbar(im)
		if detail == "some" or detail == "none":
			grid.cbar_axes[0].colorbar(im)
			grid.cbar_axes[0].toggle_label(False)
			ts = 0.001
			grid[0].text(ts, 1-ts, vtypes[0], fontsize=10, color='red', horizontalalignment='left', verticalalignment='top', transform = grid[0].transAxes)
			grid[1].text(ts, 1-ts, vtypes[1], fontsize=10, color='red', horizontalalignment='left', verticalalignment='top', transform = grid[1].transAxes)
			grid[2].text(ts, 1-ts, vtypes[2], fontsize=10, color='red', horizontalalignment='left', verticalalignment='top', transform = grid[2].transAxes)
			grid[3].text(ts, 1-ts, vtypes[3], fontsize=10, color='red', horizontalalignment='left', verticalalignment='top', transform = grid[3].transAxes)			
		if detail == "none":
			grid.cbar_axes[0].get_yaxis().set_visible(False)
		return grid












