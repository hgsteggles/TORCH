import numpy as np

class ModelData:
	fileprefix2 = "data/model_paper1/set2/"
	fileprefix3 = "data/model_paper1/set3/"
	offgrid = np.fromfile(fileprefix2 + "offgrid_times", dtype=int, sep='\n')
	injection_radii = np.fromfile("data/model_paper1/injection-sizes", dtype=float, sep='\n')
	central_densities = np.genfromtxt("data/model_paper1/central-densities")[:,1]
	central_pressures = np.genfromtxt("data/model_paper1/central-pressures")[:,1]
	size_data = np.genfromtxt("data/model_paper1/sizes", delimiter=',', skip_header=2)
	masses = [6, 9, 12, 15, 20, 30, 40, 70, 120]
	densities = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]
	times = [2.0e4, 4.0e4, 6.0e4, 8e4, 10e4]
	snapshots = [10, 20, 30, 40, 50]
	param_index = range(19, 28, 1)
	angles = [0, 30, 45, 60, 90]

	def __init__(self, useSet4=False):
		if useSet4:
			self.fileprefix3 = "data/model_paper1/set4/"

		self.mvd = np.zeros((6, 5, 9))
		self.mvt = np.zeros((6, 5, 9))

		for i in range(90):
			id = (int)(self.size_data[i,0])
			if self.size_data[i,1] == 25:
				col = (int)((id - 1) / 9)
				row = (id - 1) % 9
				self.mvd[:,col,row] = self.size_data[i,2:]
			else:
				col = (int)(self.size_data[i,1] / 10) - 1
				row = id - 19
				self.mvt[:,col,row] = self.size_data[i,2:]

		self.mvd_sw_r = self.mvd[0,:,:]
		self.mvd_sw_b = self.mvd[1,:,:]
		self.mvd_sw_t = self.mvd[2,:,:]
		self.mvd_if_r = self.mvd[3,:,:]
		self.mvd_if_b = self.mvd[4,:,:]
		self.mvd_if_t = self.mvd[5,:,:]
		self.mvt_sw_r = self.mvt[0,:,:]
		self.mvt_sw_b = self.mvt[1,:,:]
		self.mvt_sw_t = self.mvt[2,:,:]
		self.mvt_if_r = self.mvt[3,:,:]
		self.mvt_if_b = self.mvt[4,:,:]
		self.mvt_if_t = self.mvt[5,:,:]

	def getIndex(self, imass, iden):
		return len(self.masses) * iden + imass + 1

	def getAngStr(self, angle):
		if angle == 0:
			return ""
		else:
			return "_i" + str(angle)

	def getFreqStr(self, freq):
		if freq == 5:
			return ""
		elif int(freq) == freq:
			return "_f" + str(freq)
		else:
			return "_f" + str(int(10*freq)) + "e-1"

	def isSet2(self, imass, iden, itime):
		index = self.getIndex(imass, iden)
		return self.offgrid[index - 1] >= itime

	def getDataDirname(self, imass, iden, itime):
		index = self.getIndex(imass, iden)
		pad_index = "%02d" % (index,)
		isSet2 = self.isSet2(imass, iden, itime)
		filesuffix = "data_" + pad_index

		if isSet2:
			return self.fileprefix2 + filesuffix
		else:
			return self.fileprefix3 + filesuffix

	def getParamFilename(self, imass, iden, itime):
		return self.getDataDirname(imass, iden, itime) + "/parameters_" + str(self.getIndex(imass, iden)) + ".lua"

	def getRadioDirname(self, imass, iden, itime, angle, freq):
		pad_itime = "%03d" % (itime,)
		datadirname = self.getDataDirname(imass, iden, itime)
		angstr = "" if freq == 5 else self.getAngStr(angle)
		return datadirname + "/radio_" + pad_itime + angstr + self.getFreqStr(freq)

	def getDataFilename(self, imass, iden, itime):
		pad_itime = "%03d" % (itime,)
		return self.getDataDirname(imass, iden, itime) + "/data2D_" + pad_itime + ".txt"

