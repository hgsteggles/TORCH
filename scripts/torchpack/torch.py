import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
import linecache
from astropy.io import fits
import math

from mpl_toolkits.axes_grid1 import ImageGrid, Grid

def set_font_sizes(fontsize=18):
	font = {'family':'STIXGeneral','style':'normal','variant':'normal',
		'weight':'medium','size':(fontsize)}
	plt.rc('font',**font)
	plt.rc('legend',**{'fontsize':(fontsize)})
	plt.rc('xtick',**{'labelsize':fontsize})
	plt.rc('ytick',**{'labelsize':fontsize})
	plt.rc('axes',**{'labelsize':(fontsize+2)})
	#plt.rc('mathtext',**{'fontset':'custom','rm':'Bitstream Vera Sans',
	# 'it':'Bitstream Vera Sans:italic','bf':'Bitstream Vera Sans:bold'})
	plt.rc('mathtext',**{'fontset':'stix'})
	plt.rc('axes',**{'labelpad':(-0.5 + 1.5*(fontsize/5.0))})
	plt.rc('xtick.major', **{'pad':(-0.5 + 1.5*(fontsize/4.0))})
	plt.rc('ytick.major', **{'pad':(-0.5 + 1.5*(fontsize/4.0))})

set_font_sizes()

class VarType:
	name = ""
	isLog10 = False

	def __init__(self, name, units="", isLog10=False):
		self.name = name
		self.units = units
		self.isLog10 = isLog10

class Data:
	var_typenames = []

	def __init__(self, inputfile, axial=False):
		self.tsecs = float(linecache.getline(inputfile, 1))
		self.t = int(self.tsecs/31536000.0/10)*10
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
			minZ = self.get_var_raw('z').min()
			self.data = self.data[self.data[:,2] == minZ]

		if axial:
			self.mirror_data()

		self.x = []
		coords = ['x', 'y', 'z']
		for i in range(min(2, self.nd)):
			self.x.append(self.get_var_raw(coords[i]))

		self.dx = self.x[0][1] - self.x[0][0]

		self.xi = []
		self.xi.append(np.linspace(self.x[0].min(), self.x[0].max(), self.nx))
		if self.nd > 1:
			self.xi.append(np.linspace(self.x[1].min(), self.x[1].max(), self.ny))
			self.xi[0], self.xi[1] = np.meshgrid(self.xi[0], self.xi[1])

	def mirror_data(self):
		data2 = self.data[self.get_var_raw('x') != 0]
		data2[:,0] = -data2[:,0]
		self.data = np.append(self.data, data2, axis=0)

	def get_var(self, var_type):
		print "TorchData::get_var: " + var_type.name + " does not exist in this data set."
		return None

	def safe_log10(self, var):
		try:
			var[var == 0] = var[var != 0].min()
		except ValueError:
			pass
		return np.log10(var)

	def get_log10_variable(self, var_name):
		return self.safe_log10(self.get_variable(self, var_name))

	def appropriate_to_log(self, var_typename):
		return True

	def interpolate(self, var, interp):
		return scipy.interpolate.griddata((self.x[0], self.x[1]), var,
										  (self.xi[0], self.xi[1]),
										  method=interp)

class CFD_Data(Data):
	var_typenames = ['x', 'y', 'z', 'nh', 'den', 'pre', 'hii', 'nhii', 'nhi',
					 'tem', 'vel0', 'vel1', 'vel2', 'mach']

	def __init__(self, inputfile, axial=False):
		Data.__init__(self, inputfile, axial)

	def mirror_data(self):
		data2 = self.data[self.get_var_raw('x') != 0]
		data2[:,0] = -data2[:,0]
		data2[:,self.nd+3] = -data2[:,self.nd+3]
		self.data = np.append(self.data, data2, axis=0)

	def get_var(self, var_type):
		if (var_type.isLog10):
			return self.safe_log10(self.get_var_raw(var_type.name))
		else:
			return self.get_var_raw(var_type.name)
		
	def get_var_raw(self, var_typename):
		if self.nd == 0:
			return None
		elif var_typename == 'x':
			return self.data[:,0]/3.09e18
		elif var_typename == 'y' and self.nd > 1:
			return self.data[:,1]/3.09e18
		elif var_typename == 'z' and self.nd > 2:
			return self.data[:,2]/3.09e18
		elif var_typename == 'vol-cartesian':
			vol = 1
			axes = ['x', 'y', 'z']
			for i in range(self.nd):
				vol *= self.get_var_raw(axes[i])
			return vol
		elif var_typename == 'vol-cylindrical':
			rc = self.get_var_raw('x')
			return 2.0 * math.pi * rc * self.dx * self.dx
		elif var_typename == 'vol-spherical':
			rc = self.get_var_raw('x') / self.dx
			r2 = (rc + 0.5) * self.dx
			r1 = (rc - 0.5) * self.dx
			return 4.0 * math.pi * (r2 - r1) * (r2 * r2 + r1 * r2 + r1 * r1) / 3.0
		elif var_typename == 'nh':
			return self.get_var_raw('den')/1.674e-24
		elif var_typename == 'den':
			return self.data[:,self.nd]
		elif var_typename == 'pre':
			return self.data[:,self.nd+1]
		elif var_typename == 'hii':
			return self.data[:,self.nd+2]
		elif var_typename == 'nhii':
			return self.get_var_raw('hii') * self.get_var_raw('den') / 1.674e-24
		elif var_typename == 'tem':
			return (self.get_var_raw('pre') / self.get_var_raw('den')) \
				   / (8.314462e7 * (self.get_var_raw('hii') + 1))
		elif var_typename == 'ke':
			ke = 0
			for i in range(self.nd):
				ke += (self.get_var_raw('vel'+str(i)))**2
			return 0.5 * self.get_var_raw('den') * ke
		elif var_typename == 'vel0':
			return self.data[:,self.nd+3]
		elif var_typename == 'vel1' and self.nd > 1:
			return self.data[:,self.nd+4]
		elif var_typename == 'vel2' and self.nd > 2:
			return self.data[:,self.nd+5]
		elif var_typename == 'mach':
			vel = self.get_var_raw('vel0')
			mach = vel * vel
			vel = self.get_var_raw('vel1')
			if vel != None:
				mach += (vel * vel)
			vel = self.get_var_raw('vel2')
			if vel != None:
				mach += (vel * vel)
			mach = np.sqrt(mach)
			ss = np.sqrt((5.0 / 3.0) * self.get_var_raw('pre') / self.get_var_raw('den'))
			return mach / ss
		else:
			print "TorchCFD::get_var_raw: " + var_typename + \
				  " does not exist in this data set."
			return None

	def appropriate_to_log(self, var_typename):
		return var_typename != 'hii' and var_typename != 'vel0' and \
			   var_typename != 'vel1' and var_typename != 'vel2' and \
			   var_typename != 'mach'

	def convert_to_fits(self, var_type, filename):
		hdu = fits.PrimaryHDU(self.interpolate(self.get_var(var_type), "nearest"))
		hdulist = fits.HDUList([hdu])
		hdulist.writeto(filename, clobber=True)

class Cool_Data(Data):
	var_typenames = ['x', 'y', 'z', 'imlc', 'nmlc', 'recc', 'colc', 'ciec',
					 'nmc', 'fuvh', 'irh', 'crh', 'heat', 'cool', 'tot']

	def __init__(self, inputfile, axial=False):
		Data.__init__(self, inputfile, axial)

	def mirror_data(self):
		data2 = self.data[self.get_var_raw('x') != 0]
		data2[:,0] = -data2[:,0]
		self.data = np.append(self.data, data2, axis=0)

	def get_var(self, var_type):
		if (var_type.isLog10):
			return self.safe_log10(self.get_var_raw(var_type.name))
		else:
			return self.get_var_raw(var_type.name)
		
	def get_var_raw(self, var_typename):
		if self.nd != 0:
			if var_typename == 'x':
				return self.data[:,0]/3.09e18
			elif var_typename == 'y' and self.nd > 1:
				return self.data[:,1]/3.09e18
			elif var_typename == 'z' and self.nd > 2:
				return self.data[:,2]/3.09e18
			elif var_typename == 'imlc':
				return self.data[:,self.nd]
			elif var_typename == 'nmlc':
				return self.data[:,self.nd+1]
			elif var_typename == 'recc':
				return -self.data[:,self.nd+2]
			elif var_typename == 'colc':
				return self.data[:,self.nd+3]
			elif var_typename == 'ciec':
				return self.data[:,self.nd+4]
			elif var_typename == 'nmc':
				return self.data[:,self.nd+5]
			elif var_typename == 'fuvh':
				return self.data[:,self.nd+6]
			elif var_typename == 'irh':
				return self.data[:,self.nd+7]
			elif var_typename == 'crh':
				return self.data[:,self.nd+8]
			elif var_typename == 'heat':
				return self.get_var_raw('fuvh') + self.get_var_raw('irh') \
					   + self.get_var_raw('crh')
			elif var_typename == 'cool':
				return self.gat_var_raw('imlc') + self.gat_var_raw('nmlc') \
					   + self.gat_var_raw('recc') + self.gat_var_raw('colc') \
					   + self.gat_var_raw('ciec') + self.gat_var_raw('nmc')
			elif var_typename == 'lcool':
				return self.gat_var_raw('imlc') + self.gat_var_raw('recc') \
					   + self.gat_var_raw('colc') + self.gat_var_raw('ciec') \
					   + self.gat_var_raw('nmc')
			elif var_typename == 'tot':
				return self.data[:,self.nd+6] + self.data[:,self.nd+7] \
					   + self.data[:,self.nd+8] + self.data[:,self.nd] \
					   + self.data[:,self.nd+1] + self.data[:,self.nd+2] \
					   + self.data[:,self.nd+3] + self.data[:,self.nd+4] \
					   + self.data[:,self.nd+5]
		else:
			print "TorchCool::get_var: " + var_typename + \
				  " does not exist in this data set."
			return None

class InterpGrid:
	def __init__(self, x, y, nx, ny):
		self.nx = nx
		self.ny = ny
		self.x = x
		self.y = y
		self.xi = np.linspace(x.min(), x.max(), nx)
		self.yi = np.linspace(y.min(), y.max(), ny)
		self.xi, self.yi = np.meshgrid(self.xi, self.yi)

	def interpolate(self, var, interp):
		return scipy.interpolate.griddata((self.x, self.y), var,
										  (self.xi, self.yi),
										  method=interp)

class DataSetMinMax:
	def __init__(self, inputfile):
		self.min = dict()
		self.max = dict()

		minmaxfile = open(inputfile, 'r')
		for line in minmaxfile.readlines():
			columns = line.strip().split()
			if len(columns) == 6:
				self.min[columns[0]] = float(columns[3])
				self.max[columns[0]] = float(columns[4])
		minmaxfile.close()

	def get(self, key):
		return [self.min.get(key), self.max.get(key)]

class PlotParams:
	def __init__(self, datacubes, var_types, var_minmaxes, axial, interp,
				 nrows_ncols, color_maps, tight, cbar_on=True, detail="all"):
		self.datacubes = datacubes
		self.var_types = var_types
		self.var_minmaxes = var_minmaxes
		self.nrows_ncols = nrows_ncols
		self.color_maps = color_maps
		self.tight = tight
		self.cbar_on = cbar_on
		self.detail = detail
		self.axial = axial
		self.interp = interp
		self.xminmax = [None, None]
		self.yminmax = [None, None]

def imageCrop(x, y, xminmax, yminmax, var, var_min, var_max, interp, color_map, ax):
	lerp = "none"
	if interp == "linear":
		lerp = "bilinear"
	elif interp == "cubic":
		lerp = "bicubic"

	xmin = xminmax[0] if not xminmax[0] == None else x.min()
	xmax = xminmax[1] if not xminmax[1] == None else x.max()
	ymin = yminmax[0] if not yminmax[0] == None else y.min()
	ymax = yminmax[1] if not yminmax[1] == None else y.max()

	im = ax.imshow(var, vmin=var_min, vmax=var_max, origin='lower',
				   extent=[xmin,xmax,ymin,ymax],
				   interpolation=lerp, cmap=color_map)
	return im

def image(x, y, var, var_min, var_max, interp, color_map, ax):
	xminmax = [None, None]
	yminmax = [None, None]

	return imageCrop(x, y, xminmax, yminmax, var, var_min, var_max, interp, color_map, ax)


class Plotter:
	def __init__(self, nx, ny, plot_size=5, figformat='jpg', dpi=300, ncols=1):
		self.plot_size = plot_size
		self.figformat = figformat
		self.dpi = dpi
		self.cells_per_inch = nx/plot_size
		self.aspect_ratio = ny/float(nx)
		self.fig = plt.figure()

		self.ticklength = 8 * (plot_size / (5.0 * ncols))
		self.minorticklength = 4 * (plot_size / (5.0 * ncols))
		self.tickwidth = 1 * (plot_size / (5.0 * ncols))
		self.minortickwidth = 0.8 * (plot_size / (5.0 * ncols))
		self.linewidth = 2 * (plot_size / (5.0 * ncols))

	def addLog10string(self, string):
		return "\\log_{10}(" + string + ")";

	def romanise(self, string):
		return "\mathrm{" + string + "}"

	def format_units(self, var_typename):
		result = ""
		if var_typename == 'nh':
			result = "\\mathrm{cm^{-3}}"
		elif var_typename == 'den':
			result = "\\mathrm{g\\, cm^{-3}}"
		elif var_typename == 'pre':
			result = "\mathrm{Ba}"
		elif var_typename == 'nhii':
			result = "\\mathrm{cm^{-3}}"
		elif var_typename == 'tem':
			result = "\\mathrm{K}"
		elif var_typename == 'vel0' or var_typename == 'vel1' or var_typename == 'vel2':
			result = "\\mathrm{cm\\, s^{-1}}"

		if result != "":
			result = "\ [" + result + "]"

		return result

	def format_typename(self, var_typename):
		result = ""
		if var_typename == 'nh':
			result = "n_\\mathrm{H}"
		elif var_typename == 'den':
			result = "\\rho"
		elif var_typename == 'pre':
			result = "P"
		elif var_typename == 'hii':
			result = "f"
		elif var_typename == 'nhii':
			result = "n_\mathrm{HII}"
		elif var_typename == 'tem':
			result = "T"
		elif var_typename == 'vel0':
			result = "v_0"
		elif var_typename == 'vel1':
			result = "v_1"
		elif var_typename == 'vel2':
			result = "v_2"
		else:
			result = var_typename

		return result

	def format_label(self, var_type):
		var_text = self.format_typename(var_type.name)
		if var_type.isLog10:
			var_text = self.addLog10string(var_text)

		if var_type.units == "":
			var_text = var_text + self.format_units(var_type.name)
		else:
			var_text = var_text + "\ [\\mathrm{" + var_type.units + "}]"

		return "$" + var_text + "$"


	def add_quiver(self, ax, tdat, every, vmin=3.0e4, vmiddle=3.0e6):
		e = every
		ui = tdat.interpolate(tdat.get_var_raw('vel0'), 'nearest')
		vi = tdat.interpolate(tdat.get_var_raw('vel1'), 'nearest')
		#max_vel = ui.max()
		#max_vel = max(max_vel, vi.max())

		vtot = np.sqrt(ui * ui + vi * vi)
		max_vel = vtot.max()

		print max_vel

		m_high = np.logical_or(vtot < vmiddle, vtot < vmin)
		ui_high = np.ma.masked_array(ui, mask=m_high)
		vi_high = np.ma.masked_array(vi, mask=m_high)

		ax.quiver(tdat.xi[0][::e,::e], tdat.xi[1][::e,::e], ui_high[::e,::e],
				  vi_high[::e,::e], pivot='mid', units='inches', color='r',
				  scale=1.1 * max_vel * self.cells_per_inch / e)

		if vmiddle > vmin:
			m_low = np.logical_or(vtot < vmin, vtot >= vmiddle)
			ui_low = np.ma.masked_array(ui, mask=m_low)
			vi_low = np.ma.masked_array(vi, mask=m_low)

			ax.quiver(tdat.xi[0][::e,::e], tdat.xi[1][::e,::e], ui_low[::e,::e],
					  vi_low[::e,::e], pivot='mid', units='inches', color='b',
					  scale=1.1 * vmiddle * self.cells_per_inch / e)
	
	def save_plot(self, outputfile, ax=None, tight=False):
		if tight:
			self.fig.tight_layout()
		if ax == None:
			self.fig.savefig(outputfile, format=self.figformat, dpi=self.dpi,
							 bbox_inches='tight', pad_inches=0.1)
		else:
			extent = ax.window_extent().transformed(self.fig.dpi_scale_trans.inverted())
			self.fig.savefig(outputfile, format=self.figformat, dpi=self.dpi,
							 bbox_inches=extent, pad_inches=0)

	def areVarsShared(self, vtypes, vminmaxes):
		vtype = vtypes[0]
		for i in range(1, len(vtypes)):
			if vtypes[i].name != vtype.name or vtypes[i].isLog10 != vtype.isLog10:
				return False
		vmin = vminmaxes[0][0]
		vmax = vminmaxes[0][1]
		for i in range (0, len(vminmaxes)):
			if vminmaxes[i][0] == None or vminmaxes[i][0] != vmin:
				return False
			if vminmaxes[i][1] == None or vminmaxes[i][1] != vmax:
				return False
		return True

	def modifyGrid(self, grid, hasColorBar):
		for i in range(grid.ngrids):
			for axis in ['top','bottom','left','right']:
				grid[i].spines[axis].set_linewidth(self.linewidth)
				if hasColorBar:
					grid.cbar_axes[i].spines[axis].set_linewidth(self.linewidth)
			grid[i].xaxis.set_tick_params(which='major', length=self.ticklength,
										  width=self.tickwidth)
			grid[i].yaxis.set_tick_params(which='major', length=self.ticklength,
										  width=self.tickwidth)
			grid[i].xaxis.set_tick_params(which='minor', length=self.minorticklength,
										  width=self.minortickwidth)
			grid[i].yaxis.set_tick_params(which='minor', length=self.minorticklength,
										  width=self.minortickwidth)
			if hasColorBar:
				caxes = grid.cbar_axes[i]
				caxes.get_xaxis().set_tick_params(which='major', length=self.ticklength,
												  width=self.tickwidth)
				caxes.get_yaxis().set_tick_params(which='major', length=self.ticklength,
												  width=self.tickwidth)
				caxes.get_xaxis().set_tick_params(which='minor', length=self.minorticklength,
												  width=self.minortickwidth)
				caxes.get_yaxis().set_tick_params(which='minor', length=self.minorticklength,
												  width=self.minortickwidth)

	def axes1D(self, nr_nc, aspect_ratio = 1, axes_pad = 0):
		self.fig.set_size_inches(self.plot_size*nr_nc[1],
								 self.plot_size*nr_nc[0] * aspect_ratio, forward=True)
		grid = Grid(self.fig, 111, nrows_ncols = nr_nc, axes_pad = axes_pad,
					label_mode = "L", share_all = False)

		self.ticklength *= 0.8
		self.tickwidth *= 0.8
		self.minorticklength *= 0.8
		self.minortickwidth *= 0.8
		self.linewidth *= 0.8

		self.modifyGrid(grid, False)
		self.fig.set_size_inches(self.plot_size*nr_nc[1], self.plot_size*nr_nc[0] * aspect_ratio, forward=True)

		return grid

	def getGrid(self, params):
		self.fig.set_size_inches(self.plot_size*params.nrows_ncols[1],
								 self.plot_size*params.nrows_ncols[0], forward=True)
		shared_vars = False if params.var_types == None else self.areVarsShared(params.var_types, params.var_minmaxes)
		less_detail = shared_vars or params.detail == "some" or params.detail == "none"

		if less_detail:
			cbsize = str(params.nrows_ncols[1]*3.0/2.0) + "0%"
			cbpad = "0%" if params.tight else "1%"
			cbloc = "right"
			cbmode = "single" if params.cbar_on else "none"
			axpad = 0 if params.tight else 0.1*(self.plot_size / 5.0)
		else:
			cbsize = "6%"
			cbpad = "2%"
			cbloc = "right"
			cbmode = "each" if params.cbar_on else "none"
			axpad = (self.plot_size / params.nrows_ncols[1] / 5.0) * np.array([1.6, 0.4])

		return ImageGrid(self.fig, 111, nrows_ncols = params.nrows_ncols,
			axes_pad = axpad,
			label_mode = "L",
			share_all = False,
			cbar_location=cbloc,
			cbar_mode=cbmode,
			cbar_size=cbsize,
			cbar_pad=cbpad)

	def multi(self, params):
		shared_vars = self.areVarsShared(params.var_types, params.var_minmaxes)
		less_detail = shared_vars or params.detail == "some" or params.detail == "none"

		grid = self.getGrid(params)

		im = []

		for i in range(len(params.datacubes)):
			if params.detail != "none":
				if params.axial:
					grid[i].set_xlabel('$r\ \\left[\mathrm{pc}\\right]$')
					grid[i].set_ylabel('$z\ \\left[\mathrm{pc}\\right]$')
				else:
					grid[i].set_xlabel('$x\ \\left[\mathrm{pc}\\right]$')
					grid[i].set_ylabel('$y\ \\left[\mathrm{pc}\\right]$')
			else:
				grid[i].get_xaxis().set_visible(False)
				grid[i].get_yaxis().set_visible(False)

			vs = params.datacubes[i].get_var(params.var_types[i])
			vsi = params.datacubes[i].interpolate(vs, params.interp)
			vsmin = params.var_minmaxes[i][0]
			vsmax = params.var_minmaxes[i][1]

			if vsmin == None:
				vsmin = vs.min()
			if vsmax == None:
				vsmax = vs.max()

			im.append(imageCrop(params.datacubes[i].x[0], params.datacubes[i].x[1],
							params.xminmax, params.yminmax,
							vsi, vsmin, vsmax, params.interp, params.color_maps[i],
							grid[i]))

			if not less_detail:
				grid.cbar_axes[i].colorbar(im[i])
				cbax = grid.cbar_axes[i]
				cblax = cbax.axis[cbax.orientation]
				cblax.label.set_text(self.format_label(params.var_types[i]))

		if less_detail:
			grid.cbar_axes[0].colorbar(im[0])
			cbax = grid.cbar_axes[0]
			cblax = cbax.axis[cbax.orientation]
			cblax.label.set_text(self.format_label(params.var_types[0]))

		if params.detail == "some" or params.detail == "none":
			grid.cbar_axes[0].colorbar(im[0])
			grid.cbar_axes[0].toggle_label(False)
			ts = 0.04

			for i in range(len(params.var_types)):
				grid[i].text(ts, 1-ts, self.format_label(params.var_types[i]),
							 color='red', horizontalalignment='left',
							 verticalalignment='top', transform=grid[i].transAxes)

		if params.detail == "none":
			grid.cbar_axes[0].get_yaxis().set_visible(False)

		self.modifyGrid(grid, True)

		return grid

	def multi2(self, params, grid):
		shared_vars = self.areVarsShared(params.var_types, params.var_minmaxes)
		less_detail = shared_vars or params.detail == "some" or params.detail == "none"

		im = []

		for i in range(len(params.datacubes)):
			if params.detail != "none":
				if params.axial:
					grid[i].set_xlabel('$r\ \\left[\mathrm{pc}\\right]$')
					grid[i].set_ylabel('$z\ \\left[\mathrm{pc}\\right]$')
				else:
					grid[i].set_xlabel('$x\ \\left[\mathrm{pc}\\right]$')
					grid[i].set_ylabel('$y\ \\left[\mathrm{pc}\\right]$')
			else:
				grid[i].get_xaxis().set_visible(False)
				grid[i].get_yaxis().set_visible(False)

			vs = params.datacubes[i].get_var(params.var_types[i])
			vsi = params.datacubes[i].interpolate(vs, params.interp)
			vsmin = params.var_minmaxes[i][0]
			vsmax = params.var_minmaxes[i][1]

			if vsmin == None:
				vsmin = vs.min()
			if vsmax == None:
				vsmax = vs.max()

			im.append(imageCrop(params.datacubes[i].x[0], params.datacubes[i].x[1],
							params.xminmax, params.yminmax,
							vsi, vsmin, vsmax, params.interp, params.color_maps[i],
							grid[i]))

			xmin = min(params.datacubes[i].x[0])
			xmax = max(params.datacubes[i].x[0])
			ymin = min(params.datacubes[i].x[1])
			ymax = max(params.datacubes[i].x[1])

			grid[i].set_xlim([xmin, xmax])
			grid[i].set_ylim([ymin, ymax])

			if not less_detail:
				cbax = grid.cbar_axes[i]
				cb = grid.fig.colorbar(im[i], cax=cbax, orientation='horizontal')
				cbax.xaxis.label.set_text(self.format_label(params.var_types[i]))

		if less_detail:
			grid.cbar_axes[0].colorbar(im[0])
			cbax = grid.cbar_axes[0]
			cblax = cbax.axis[cbax.orientation]
			cblax.label.set_text(self.format_label(params.var_types[0]))

		if params.detail == "some" or params.detail == "none":
			grid.cbar_axes[0].colorbar(im[0])
			grid.cbar_axes[0].toggle_label(False)
			ts = 0.04

			for i in range(len(params.var_types)):
				grid[i].text(ts, 1-ts, self.format_label(params.var_types[i]),
							 color='red', horizontalalignment='left',
							 verticalalignment='top', transform = grid[i].transAxes)

		if params.detail == "none":
			grid.cbar_axes[0].get_yaxis().set_visible(False)

		self.modifyGrid(grid, True)

		return grid

	def histstep(self, ax, data, bins, errorcentres=None, normed=False, **kwargs):
		#grid[0].hist(data, bins=bins)
		histdata, bin_edges = np.histogram(data, bins=bins)
		N = sum(histdata)

		hist = histdata if not normed else histdata / float(N)

		ax.step(bin_edges, np.concatenate((hist, [0]), axis=0), where='post', **kwargs)

		if errorcentres is not None:
			del kwargs['label']
			yerr = np.sqrt(histdata)
			if normed:
				yerr = yerr / float(N)

			ax.errorbar(errorcentres, hist, yerr=yerr, fmt='.', **kwargs)

class FancyAxesGrid:
	def __init__(self, width_ratios, height_ratios, hspace = 0.01, vspace = 0.01,
				 cspace = 0.01, csize = 0.01, cpad=0.03, fig_width=5,
				 border=(0.01, 0.01, 0.01, 0.01), fig_format='png', dpi=300):
		self.ncols_nrows = (len(width_ratios), len(height_ratios))
		self.hspace = hspace
		self.vspace = vspace
		self.cspace = cspace
		self.csize = csize
		self.cpad = cpad
		self.border = border
		self.fig_format = fig_format
		self.dpi = dpi
		self.fig = plt.figure()
		self.ticklength = 2
		self.visible = []
		self.ngrids = len(width_ratios) * len(height_ratios)

		# Initialise all visibles
		for i in range(self.ncols_nrows[0]):
			vis = []
			for j in range(self.ncols_nrows[1]):
				vis.append(True)
			self.visible.append(vis)

		# Calculate fig, plot, and axes sizes and fractions
		self.ax_wfractions = []
		self.ax_hfractions = []

		plot_w = fig_width*(1.0 - border[0] - border[1])
		plot_frac_x = (plot_w)/float(fig_width)
		denom_w = sum(width_ratios)/(1.0 - (self.ncols_nrows[0] - 1) * hspace)

		for i in range(len(width_ratios)):
			self.ax_wfractions.append(width_ratios[i] * plot_frac_x / denom_w)

		sum_ifig_w = 0

		for j in range(len(height_ratios)):
			ifig_rat = height_ratios[j]/float(width_ratios[0])
			sum_ifig_w = sum_ifig_w + ifig_rat * self.ax_wfractions[0] * plot_w

		plot_h = sum_ifig_w / (1.0 - self.ncols_nrows[1] * (cspace + csize + vspace))

		self.fig_size = (fig_width, plot_h / (1.0 - border[2] - border[3]))
		self.fig.set_size_inches(self.fig_size[0], self.fig_size[1])
		self.plot_size = (plot_w, plot_h)

		plot_frac_y = self.plot_size[1]/float(self.fig_size[1])
		denom_h = sum(height_ratios)/(1.0 - self.ncols_nrows[1] * (vspace + cspace + csize))

		for i in range(len(height_ratios)):
			self.ax_hfractions.append(height_ratios[i] * plot_frac_y / denom_h)


		# Set up axes
		posX = self.border[0]

		grid = []
		cgrid = []

		for col in range(self.ncols_nrows[0]):
			grid_col = []
			cgrid_col = []
			posY = 1.0 - self.border[3]

			for row in range(self.ncols_nrows[1]):
				posY = posY - self.ax_hfractions[row] - vspace - csize - cspace

				ax = plt.axes([posX, posY, self.ax_wfractions[col], self.ax_hfractions[row]],
							  autoscale_on=False)
				grid_col.append(ax)
				cax = plt.axes([posX + (self.cpad * self.ax_wfractions[col]),
								posY + self.ax_hfractions[row] + cspace + (self.cpad * csize),
								self.ax_wfractions[col] * (1.0 - (2.0 * self.cpad)),
								csize * (1.0 - (2.0 * self.cpad))],
							   autoscale_on=False, axisbelow=False)
				cgrid_col.append(cax)

			posX = posX + self.ax_wfractions[col] + hspace

			grid.append(grid_col)
			cgrid.append(cgrid_col)

		self.grid = tuple(grid)
		self.cgrid = tuple(cgrid)

		self.cbar_axes = []
		for row in range(self.ncols_nrows[1]):
			for col in range(self.ncols_nrows[0]):
				self.cbar_axes.append(self.cgrid[col][row])

	def __getitem__(self, item):
		return self.grid[item % self.ncols_nrows[0]][int(item / self.ncols_nrows[0])]

	def set_visible(self, i, j, isVisible):
		self.visible[i][j] = isVisible
		if isVisible:
			self.cgrid[i][j].set_axis_on()
			self.grid[i][j].set_axis_on()
		else:
			self.cgrid[i][j].set_axis_off()
			self.grid[i][j].set_axis_off()

	def is_visible(self, i, j):
		return self.visible[i][j]

	def remove_ticklabels(self):
		# Remove tick labels.
		for col in range(self.ncols_nrows[0]):
			for row in range(self.ncols_nrows[1]):
				if col != 0 and self.visible[col-1][row]:
					plt.setp(self.grid[col][row].get_yticklabels(), visible=False)
					self.grid[col][row].set_ylabel("")
				if row != self.ncols_nrows[1] - 1 and self.visible[col][row+1]:
					plt.setp(self.grid[col][row].get_xticklabels(), visible=False)
					self.grid[col][row].set_xlabel("")

	def update_cbar(self):
		# Move cbar xaxis to top.
		for row in range(self.ncols_nrows[1]):
			for col in range(self.ncols_nrows[0]):
				self.cgrid[col][row].get_xaxis().set_ticks_position('top')
				self.cgrid[col][row].xaxis.set_tick_params(length=self.ticklength)
				self.cgrid[col][row].xaxis.set_label_position('top')

	def update(self):
		self.remove_ticklabels()

		self.update_cbar()

	def save_plot(self, outputfile):
		self.fig.savefig(outputfile, format=self.fig_format, dpi=self.dpi, pad_inches=0)

