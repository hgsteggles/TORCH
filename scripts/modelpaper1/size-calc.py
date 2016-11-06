import numpy as np
import math
import sys
import lupa
from scipy import stats

lua = lupa.LuaRuntime(unpack_returned_tuples=True)

def load_src(name, fpath):
    import os, imp
    return imp.load_source(name, os.path.join(os.path.dirname(__file__), fpath))

load_src("koo", "../torchpack/koo.py")

import koo

load_src("mp1", "../torchpack/modelpaper1.py")
import mp1
mp1_data = mp1.ModelData()

print_sizes = False
print_mvd_table = False
print_mvt_table = False

data = np.genfromtxt("data/model_paper1/refdata/zams.txt", skip_header=1)

massList = data[:,0]
logQList = data[:,1]
mdotList = data[:,2]
vinfList = data[:,3]

CM2PC = 3.24078e-19
S2YR = 3.17098e-8
YR2S = 3.154e7
mh = 1.6738232e-24
R_GAS = 8.3144598e7
KB = 1.38e-16

mvd_nhs = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]
mvt_nhs = [3.2e4, 3.2e4, 3.2e4, 3.2e4, 3.2e4]
mvt_times = [2.0e4, 4.0e4, 6.0e4, 8.0e4, 10.0e4]
mvd_times = [5.0e4, 5.0e4, 5.0e4, 5.0e4, 5.0e4]

######### MVD

#### High Temp Wind

xrange = []
xrange.append([15, 32, 62, 90, 134, 190, 260, 380, 800])
xrange.append([10, 22, 48, 80, 110, 140, 190, 300, 500])
xrange.append([6, 15, 34, 60, 90, 120, 150, 280, 400])
xrange.append([4, 9, 24, 46, 66, 100, 120, 180, 250])
xrange.append([2.52, 6, 16, 30, 54, 80, 100, 150, 190])

for i in range(5):
	for j in range(9):
		xrange[i][j] = xrange[i][j] * 0.1 / 15.0

mvd_xrange = []
mvd_xrange.append([0.10, 0.24, 0.46, 0.70, 0.90, 1.26, 1.73, 2.53, 5.33])
mvd_xrange.append([0.06, 0.15, 0.38, 0.60, 0.80, 0.90, 1.27, 2.00, 3.33])
mvd_xrange.append([0.04, 0.11, 0.25, 0.45, 0.65, 0.90, 1.00, 1.87, 2.67])
mvd_xrange.append([0.03, 0.06, 0.20, 0.35, 0.55, 0.80, 1.00, 1.2, 1.67])
mvd_xrange.append([0.02, 0.05, 0.14, 0.25, 0.45, 0.60, 0.80, 1, 1.27])
mvd_xrange = np.array(mvd_xrange, np.float64)

mvd_ref = 34.0

mvd_hsize = []
mvd_hsize.append([16.5, 4.25, 5.50, 5.00, 6.00, 8.00, 11.0, 29.0, 24.0])
mvd_hsize.append([19.5, 8.25, 4.00, 5.00, 6.25, 10.5, 13.0, 27.0, 24.0])
mvd_hsize.append([15.0, 10.5, 6.75, 6.00, 6.50, 9.50, 12.2, 23.0, 27.0])
mvd_hsize.append([16.5, 7.00, 7.50, 6.00, 7.00, 11.0, 18.5, 24.5, 25.5])
mvd_hsize.append([14.5, 5.00, 8.50, 5.50, 7.50, 16.0, 15.0, 25.0, 26.5])
mvd_hsize = np.array(mvd_hsize, np.float64)

mvd_vup = []
mvd_vup.append([6.00, 2.25, 3.00, 1.50, 1.50, 2.00, 2.50, 3.00, 7.50])
mvd_vup.append([5.50, 5.50, 1.00, 2.50, 2.00, 3.00, 3.00, 3.50, 3.50])
mvd_vup.append([5.00, 3.00, 3.00, 3.00, 3.00, 3.00, 5.00, 4.00, 4.00])
mvd_vup.append([7.00, 2.00, 3.50, 2.50, 4.00, 3.00, 4.00, 5.00, 5.00])
mvd_vup.append([6.00, 2.00, 3.00, 3.00, 3.25, 4.00, 4.00, 6.00, 6.00])
mvd_vup = np.array(mvd_vup, np.float64)

mvd_vdown = []
mvd_vdown.append([11.0, 12.0, 10.5, 10.5, 13.0, 17.5, 18.5, 15.0, 14.0])
mvd_vdown.append([10.0, 4.00, 10.0, 9.00, 11.0, 14.5, 16.5, 15.5, 15.5])
mvd_vdown.append([10.0, 7.00, 8.50, 9.00, 10.0, 13.5, 15.0, 14.0, 16.0])
mvd_vdown.append([10.5, 8.00, 6.50, 8.00, 8.50, 10.5, 13.5, 16.0, 16.5])
mvd_vdown.append([10.5, 4.50, 5.00, 9.00, 8.50, 9.50, 12.5, 16.0, 19.0])
mvd_vdown = np.array(mvd_vdown, np.float64)


### Ionisation Front

mvd_i_hsize = []
mvd_i_hsize.append([25.0, 24.0, 25.0, 25.0, 30.0, 30.0, 29.0, 32.0, 25.0])
mvd_i_hsize.append([26.0, 27.0, 22.0, 23.0, 27.0, 32.5, 30.0, 29.5, 25.5])
mvd_i_hsize.append([23.0, 24.0, 24.0, 23.0, 25.5, 26.0, 30.0, 28.0, 29.0])
mvd_i_hsize.append([23.5, 26.5, 21.0, 21.5, 23.5, 22.5, 23.0, 27.0, 29.0])
mvd_i_hsize.append([21.0, 19.5, 19.5, 21.5, 22.0, 23.5, 22.5, 26.5, 29.0])
mvd_i_hsize = np.array(mvd_i_hsize, np.float64)

mvd_i_vup = []
mvd_i_vup.append([11.5, 11.0, 11.0, 9.50, 10.0, 8.00, 6.00, 4.50, 8.50])
mvd_i_vup.append([10.5, 12.0, 8.50, 9.00, 10.0, 10.0, 8.00, 5.50, 4.00])
mvd_i_vup.append([8.50, 12.0, 10.5, 10.5, 10.0, 8.50, 9.00, 6.00, 5.00])
mvd_i_vup.append([11.0, 12.0, 10.5, 9.00, 10.5, 9.00, 8.00, 7.50, 6.50])
mvd_i_vup.append([9.50, 9.50, 8.50, 9.50, 10.0, 10.0, 8.50, 8.00, 7.50])
mvd_i_vup = np.array(mvd_i_vup, np.float64)

mvd_i_vdown = []
mvd_i_vdown.append([14.5, 15.0, 14.0, 15.0, 18.0, 21.0, 21.0, 18.5, 15.0])
mvd_i_vdown.append([14.5, 16.0, 13.5, 13.5, 16.0, 22.0, 20.0, 17.5, 16.0])
mvd_i_vdown.append([13.5, 13.0, 15.0, 13.5, 15.5, 16.0, 19.0, 16.0, 17.0])
mvd_i_vdown.append([12.0, 14.0, 12.0, 13.0, 13.5, 14.5, 16.0, 17.0, 18.0])
mvd_i_vdown.append([13.0, 10.5, 11.0, 13.0, 12.5, 13.5, 14.0, 16.0, 20.0])
mvd_i_vdown = np.array(mvd_i_vdown, np.float64)


######### MVT

### High Temp Wind

mvt_xrange = []
mvt_xrange.append([0.040, 0.10, 0.20, 0.28, 0.40, 0.52, 0.7, 1.0, 2.2])
mvt_xrange.append([0.035, 0.09, 0.22, 0.32, 0.50, 0.80, 1.0, 1.9, 2.5])
mvt_xrange.append([0.040, 0.10, 0.24, 0.40, 0.65, 0.90, 1.0, 2.2, 3.0])
mvt_xrange.append([0.040, 0.08, 0.25, 0.48, 0.74, 1.00, 1.3, 2.3, 3.2])
mvt_xrange.append([0.040, 0.10, 0.25, 0.54, 0.80, 1.00, 1.6, 2.6, 3.6])
mvt_xrange = np.array(mvt_xrange, np.float64)

mvt_ref = 34.0

mvt_hsize = []
mvt_hsize.append([19.0, 7.00, 5.50, 8.00, 7.00, 11.5, 20.0, 22.0, 11.0])
mvt_hsize.append([20.0, 6.50, 5.50, 8.00, 8.00, 8.50, 12.0, 20.0, 24.5])
mvt_hsize.append([17.5, 8.00, 6.00, 7.00, 8.50, 9.00, 14.0, 21.5, 24.0])
mvt_hsize.append([19.0, 7.00, 7.50, 5.50, 6.00, 8.00, 11.0, 19.0, 26.0])
mvt_hsize.append([19.0, 6.00, 7.50, 5.50, 5.00, 8.00, 11.0, 18.5, 26.0])
mvt_hsize = np.array(mvt_hsize, np.float64)

mvt_vup = []
mvt_vup.append([7.50, 3.00, 3.00, 5.00, 3.00, 6.00, 6.00, 6.00, 5.00])
mvt_vup.append([6.00, 4.00, 2.00, 4.00, 4.00, 3.00, 5.00, 4.00, 4.00])
mvt_vup.append([5.00, 5.00, 4.00, 3.00, 3.00, 3.00, 5.00, 3.00, 3.00])
mvt_vup.append([8.00, 4.50, 3.00, 3.00, 2.00, 3.00, 3.00, 3.00, 3.50])
mvt_vup.append([8.00, 6.00, 6.00, 3.00, 2.00, 2.50, 3.00, 3.00, 3.50])
mvt_vup = np.array(mvt_vup, np.float64)

mvt_vdown = []
mvt_vdown.append([11.0, 7.50, 8.50, 6.00, 9.00, 12.0, 13.5, 16.0, 9.00])
mvt_vdown.append([9.50, 10.0, 12.0, 10.5, 12.5, 14.0, 14.0, 12.5, 15.5])
mvt_vdown.append([11.0, 7.50, 7.50, 12.0, 12.0, 14.0, 14.0, 12.5, 15.0])
mvt_vdown.append([10.5, 3.00, 11.0, 12.5, 13.0, 14.0, 17.0, 13.0, 15.0])
mvt_vdown.append([11.0, 2.00, 9.00, 12.0, 14.0, 17.5, 15.5, 14.0, 15.0])
mvt_vdown = np.array(mvt_vdown, np.float64)


### Ionisation Front

mvt_i_hsize = []
mvt_i_hsize.append([27.0, 24.5, 23.5, 27.0, 26.0, 27.0, 25.0, 24.5, 21.5])
mvt_i_hsize.append([28.0, 31.0, 27.0, 31.5, 30.0, 25.5, 16.5, 24.0, 26.0])
mvt_i_hsize.append([25.5, 27.5, 27.5, 29.0, 27.5, 27.5, 32.0, 25.0, 25.5])
mvt_i_hsize.append([28.0, 33.0, 28.0, 27.0, 28.0, 29.0, 28.5, 29.5, 27.0])
mvt_i_hsize.append([27.5, 26.0, 29.0, 26.0, 28.0, 32.0, 26.0, 30.0, 27.5])
mvt_i_hsize = np.array(mvt_i_hsize, np.float64)

mvt_i_vup = []
mvt_i_vup.append([12.5, 13.0, 10.0, 13.0, 11.0, 11.0, 9.00, 8.00, 5.00])
mvt_i_vup.append([12.0, 15.5, 12.0, 14.0, 12.0, 9.00, 8.00, 6.00, 5.50])
mvt_i_vup.append([10.0, 14.0, 12.0, 13.0, 11.0, 9.00, 9.00, 4.50, 5.00])
mvt_i_vup.append([13.5, 16.0, 12.5, 12.0, 11.0, 9.00, 7.50, 5.00, 4.50])
mvt_i_vup.append([12.0, 13.0, 13.0, 11.5, 11.0, 10.0, 6.50, 4.50, 4.50])
mvt_i_vup = np.array(mvt_i_vup, np.float64)

mvt_i_vdown = []
mvt_i_vdown.append([16.0, 14.0, 13.5, 14.0, 14.5, 16.0, 15.0, 17.0, 12.5])
mvt_i_vdown.append([14.5, 16.5, 16.5, 17.0, 16.5, 16.0, 16.0, 12.5, 17.0])
mvt_i_vdown.append([16.5, 14.5, 16.5, 16.5, 16.0, 17.0, 20.0, 13.5, 15.0])
mvt_i_vdown.append([16.5, 16.5, 16.5, 16.0, 16.0, 17.0, 19.0, 14.0, 16.0])
mvt_i_vdown.append([15.5, 17.0, 17.0, 15.5, 17.0, 21.5, 17.5, 14.0, 15.0])
mvt_i_vdown = np.array(mvt_i_vdown, np.float64)

mvt_param_index = range(19, 28, 1)
mvt_snapshots = [10, 20, 30, 40, 50]

mvd_SW_r = np.multiply(mvd_xrange, mvd_hsize) * 0.5 / mvd_ref
mvd_SW_t = mvd_xrange * mvd_vup / mvd_ref
mvd_SW_b = mvd_xrange * mvd_vdown / mvd_ref

mvd_IF_r = np.multiply(mvd_xrange, mvd_i_hsize) * 0.5 / mvd_ref
mvd_IF_t = np.multiply(mvd_xrange, mvd_i_vup) / mvd_ref
mvd_IF_b = np.multiply(mvd_xrange, mvd_i_vdown) / mvd_ref

mvt_SW_r = np.multiply(mvt_xrange, mvt_hsize) * 0.5 / mvt_ref
mvt_SW_t = mvt_xrange * mvt_vup / mvt_ref
mvt_SW_b = mvt_xrange * mvt_vdown / mvt_ref

mvt_IF_r = np.multiply(mvt_xrange, mvt_i_hsize) * 0.5 / mvt_ref
mvt_IF_t = np.multiply(mvt_xrange, mvt_i_vup) / mvt_ref
mvt_IF_b = np.multiply(mvt_xrange, mvt_i_vdown) / mvt_ref

### Injection sizes

snapshots = [10, 20, 30, 40, 50]

def getInjectionSizesMVD():
	sizes = np.zeros((5, 9))

	for i in range(45):
		col = int(i / 9)
		row = i % 9

		filename = mp1_data.getParamFilename(row, col, 25)
		filestring = open(filename, 'r').read()

		table = lua.eval("{" + filestring + "}")
		if table == None:
			print filestring
			sys.exit()
		ncellsx = table["Parameters"]["Grid"]["no_cells_x"]
		length = table["Parameters"]["Grid"]["side_length"]

		sizes[col][row] = (10.0 / ncellsx) * length * CM2PC

	return sizes

def getInjectionSizesMVT():
	sizes = np.zeros((5, 9))

	for col in range(5):
		for row in range(9):
			filename = mp1_data.getParamFilename(row, 2, snapshots[col])
			filestring = open(filename, 'r').read()

			table = lua.eval("{" + filestring + "}")
			if table == None:
				print filestring
				sys.exit()
			ncellsx = table["Parameters"]["Grid"]["no_cells_x"]
			length = table["Parameters"]["Grid"]["side_length"]

			sizes[col][row] = (10.0 / ncellsx) * length * CM2PC

	return sizes

mvd_inj = getInjectionSizesMVD()
mvt_inj = getInjectionSizesMVT()

### Ionized ambient gas density

#mvd_iden = []
#mvd_iden.append([130.0,   150.0,  250.0,  400.0,  650.0,  900.0, 1200.0, 1800.0,  2200.0])
#mvd_iden.append([280.0,   250.0,  380.0,  560.0, 1000.0, 1300.0, 1800.0, 2000.0,  2200.0])
#mvd_iden.append([610.0,   480.0,  650.0,  900.0, 1300.0, 1800.0, 2200.0, 2400.0,  2700.0])
#mvd_iden.append([800.0,   900.0, 1000.0, 1400.0, 2100.0, 2500.0, 3000.0, 5000.0,  5000.0])
#mvd_iden.append([2000.0, 2200.0, 2200.0, 2400.0, 3000.0, 4500.0, 5500.0, 7000.0, 12000.0])

#mvt_iden = []
#mvt_iden.append([520.0, 620.0, 1000.0, 1700.0, 2600.0, 4500.0, 6500.0, 9000.0, 12000.0])
#mvt_iden.append([620.0, 520.0,  750.0, 1000.0, 1600.0, 2400.0, 2500.0, 5000.0,  8000.0])
#mvt_iden.append([580.0, 520.0,  600.0,  800.0, 1400.0, 1800.0, 1900.0, 2500.0,  3000.0])
#mvt_iden.append([520.0, 550.0,  550.0,  700.0, 1100.0, 1300.0, 1500.0, 1700.0,  1800.0])
#mvt_iden.append([520.0, 540.0,  520.0,  620.0,  850.0, 1200.0, 1300.0, 1400.0,  1500.0])

mvd_iden = []
for i in range(5):
	arr = []
	for j in range(9):
		R_strom = koo.calcStromgrenRadius(logQList[j], mvd_nhs[i]) * CM2PC
		#R_strom = max(R_strom, mvd_inj[i][j])
		arr.append(mvd_nhs[i] * math.pow(R_strom / mvd_IF_r[i][j], 1.5))
	mvd_iden.append(arr)

mvt_iden = []
for i in range(5):
	arr = []
	for j in range(9):
		R_strom = koo.calcStromgrenRadius(logQList[j], mvt_nhs[i]) * CM2PC
		#R_strom = max(R_strom, mvt_inj[i][j])
		arr.append(mvt_nhs[i] * math.pow(R_strom / mvt_IF_r[i][j], 1.5))
	mvt_iden.append(arr)


if print_sizes:
	for index in range(45):
		col = int(index / 9)
		row = index % 9

		SW_r = "{:6.5f}".format(mvd_SW_r[col][row])
		SW_b = "{:6.5f}".format(mvd_SW_b[col][row])
		SW_t = "{:6.5f}".format(mvd_SW_t[col][row])
		IF_r = "{:6.5f}".format(mvd_IF_r[col][row])
		IF_b = "{:6.5f}".format(mvd_IF_b[col][row])
		IF_t = "{:6.5f}".format(mvd_IF_t[col][row])

		print str(index + 1) + ",25," +  \
			  SW_r + "," + SW_b + "," + SW_t + "," + \
			  IF_r + "," + IF_b + "," + IF_t

	for col in range(5):
		for row in range(9):
			SW_r = "{:6.5f}".format(mvt_SW_r[col][row])
			SW_b = "{:6.5f}".format(mvt_SW_b[col][row])
			SW_t = "{:6.5f}".format(mvt_SW_t[col][row])
			IF_r = "{:6.5f}".format(mvt_IF_r[col][row])
			IF_b = "{:6.5f}".format(mvt_IF_b[col][row])
			IF_t = "{:6.5f}".format(mvt_IF_t[col][row])
			print str(mvt_param_index[row]) + "," + str(mvt_snapshots[col]) + "," +  \
				  SW_r + "," + SW_b + "," + SW_t + "," + \
				  IF_r + "," + IF_b + "," + IF_t

def windHiiVelRatio(isMVD):
	for row in range(9):
		mdot = mdotList[row]
		vinf = vinfList[row]
		times = mvd_times if isMVD else mvt_times

		tableline =  str(massList[row])
		for col in range(5):
			nh = mvd_iden[col][row] if isMVD else  mvt_iden[col][row]
			#nh = mvd_nhs[col] if isMVD else  mvt_nhs[col]

			R = koo.calcBubbleRadiusRB(mdot, vinf, mh * nh, times[col] * YR2S)
			v = R / (times[col] * YR2S)
			ci = koo.soundSpeed(8000.0, 0.5)
			tableline = tableline + " & " + str(v / ci)

		print tableline + "\\\\"

def printCritVel(isMVD):
	for row in range(9):
		mdot = mdotList[row]
		vinf = vinfList[row]

		tableline =  str(massList[row])
		for col in range(5):
			nh = mvd_iden[col][row] if isMVD else  mvt_iden[col][row]
			nh = mvd_nhs[col] if isMVD else  mvt_nhs[col]

			tableline = tableline + " & "
			R = vinf / koo.calcCritVel(mdot, vinf, mh * nh)
			tableline = tableline + str(R)
		print tableline + " \\\\"

def printAnalyticalTimes(isMVD):
	T = 8000.0
	mu = 0.5

	for row in range(9):
		mass = massList[row]
		mdot = mdotList[row]
		vinf = vinfList[row]
		nh = mvd_iden if isMVD else  mvt_iden
		times = mvd_times if isMVD else mvt_times

		tableline = "\multirow{3}{*}{\\num{" + str(mass) + "}} & {$t_\\mathrm{P}^\\mathrm{RB}$}"
		for col in range(5):
			tableline = tableline + " & "
			R = koo.calcConfineTimeRB(mdot, vinf, mh*nh[col][row], T, mu) * S2YR / 1000.0
			tableline = tableline + str(R)
		print tableline + " \\\\"

		tableline = "& {$t_\\mathrm{P}^\\mathrm{PRB}$}"
		for col in range(5):
			tableline = tableline + " & "
			R = koo.calcConfineTimePRB(mdot, vinf, mh*nh[col][row], T, mu) * S2YR / 1000.0
			tableline = tableline + str(R)
		print tableline + " \\\\"

		tableline = "& {$t_\\mathrm{P}^\\mathrm{AB}$}"
		for col in range(5):
			tableline = tableline + " & "
			R = koo.calcConfineTimeAB(mdot, vinf, mh*nh[col][row], T, mu) * S2YR / 1000.0
			tableline = tableline + str(R)
		print tableline + " \\\\"

		if row != 8:
			print "& & & & & & \\\\"

def printAnalyticalSizes(r, isMVD):
	T = 300.0
	mu = 1.0

	for row in range(9):
		mass = massList[row]
		mdot = mdotList[row]
		vinf = vinfList[row]
		nh = mvd_iden if isMVD else  mvt_iden
		times = mvd_times if isMVD else mvt_times

		tableline =  "\multirow{5}{*}{\\num{" + str(mass) + "}} & {$D_\mathrm{s}(\\frac{\\pi}{2})$}"
		for col in range(5):
			tableline = tableline + " & "
			tableline = tableline + str(r[col][row])
		print tableline + " \\\\"

		tableline = "& {$R_\\mathrm{s}^\\mathrm{RB}$}"
		for col in range(5):
			tableline = tableline + " & "
			R = koo.calcFinalRadiusRB(mdot, vinf, mh*nh[col][row], T, mu) * CM2PC
			tableline = tableline + str(R)
		print tableline + " \\\\"

		tableline = "& {$R_\\mathrm{s}^\\mathrm{PRB}$}"
		for col in range(5):
			tableline = tableline + " & "
			R = koo.calcFinalRadiusPRB(mdot, vinf, mh*nh[col][row], T, mu) * CM2PC
			tableline = tableline + str(R)
		print tableline + " \\\\"

		tableline = "& {$R_\\mathrm{s}^\\mathrm{FAB}$}"
		for col in range(5):
			tableline = tableline + " & "
			R = koo.calcFinalRadiusFAB(mdot, vinf, mh*nh[col][row], T, mu) * CM2PC
			tableline = tableline + str(R)
		print tableline + " \\\\"

		tableline = "& {$R_\\mathrm{s}^\\mathrm{PAB}$}"
		for col in range(5):
			tableline = tableline + " & "
			R = koo.calcFinalRadiusPAB(mdot, vinf, mh*nh[col][row], T, mu) * CM2PC
			tableline = tableline + str(R)
		print tableline + " \\\\"

		if row != 8:
			print "& & & & & & \\\\"

def printAdiabaticTable(r, isMVD):
	T = 300.0
	mu = 1.0

	diststr = "D_\\mathrm{s}"
	numstr = "2"

	for row in range(9):
		mass = massList[row]
		mdot = mdotList[row]
		vinf = vinfList[row]
		logQ = logQList[row]
		nhs_hii = mvd_iden if isMVD else  mvt_iden
		times = mvd_times if isMVD else mvt_times
		nhs_init = mvd_nhs if isMVD else mvt_nhs
		inj = mvd_inj if isMVD else mvt_inj

		tableline = "\multirow{" + numstr + "}{*}{\\num{" + str(mass) + "}} & {$" + diststr + "(\\frac{\\pi}{2})$}"
		for col in range(5):
			tableline = tableline + " & "
			tableline = tableline + str(r[col][row])
		print tableline + " \\\\"

		tableline = "& {$R_\\mathrm{P}^\\mathrm{PAB}$}"
		for col in range(5):
			tableline = tableline + " & "
			tinj = koo.calcInjectionTimeAB(mdot, vinf, mh*nhs_init[col], inj[col][row] / CM2PC)
			tinj = 0
			#print tinj * S2YR
			t = (koo.calcConfineTimeAB(mdot, vinf, mh*nhs_init[col], T, mu) - tinj) * S2YR
			if times[col] <= t:
				R = koo.calcBubbleRadiusFAB(mdot, vinf, mh*nhs_init[col], times[col] * YR2S + tinj) * CM2PC
			else:
				R = koo.calcFinalRadiusFAB(mdot, vinf, mh*nhs_init[col], T, mu) * CM2PC
			tableline = tableline + str(R)
		print tableline + " \\\\"

		if row != 8:
			print "& & & & & & \\\\"

def printPaperTable(r, t, b, isWind, isMVD):
	T = 8000.0
	mu = 0.5

	diststr = "D_\\mathrm{s}" if isWind else "D_\\mathrm{i}"
	numstr = "4" if isWind else "4"

	for row in range(9):
		mass = massList[row]
		mdot = mdotList[row]
		vinf = vinfList[row]
		logQ = logQList[row]
		nhs_hii = mvd_iden if isMVD else  mvt_iden
		times = mvd_times if isMVD else mvt_times
		nhs_init = mvd_nhs if isMVD else mvt_nhs
		inj = mvd_inj if isMVD else mvt_inj

		tableline = "\multirow{" + numstr + "}{*}{\\num{" + str(mass) + "}} & {$" + diststr + "(0)$}"
		for col in range(5):
			tableline = tableline + " & "
			tableline = tableline + str(t[col][row])
		print tableline + " \\\\"

		tableline =  "& {$" + diststr + "(\\frac{\\pi}{2})$}"
		for col in range(5):
			tableline = tableline + " & "
			tableline = tableline + str(r[col][row])
		print tableline + " \\\\"

		tableline = "& {$" + diststr + "(\\pi)$}"
		for col in range(5):
			tableline = tableline + " & "
			tableline = tableline + str(b[col][row])
		print tableline + " \\\\"

		if isWind:
			tableline = "& {$R_\\mathrm{P}$}"
			for col in range(5):
				tableline = tableline + " & "
				R = koo.calcFinalRadiusRB(mdot, vinf, mh*nhs_hii[col][row], T, mu) * CM2PC
				tableline = tableline + str(R)
			print tableline + " \\\\"
		else:
			tableline = "& {$R_\\mathrm{raga-1}$}"
			for col in range(5):
				tableline = tableline + " & "
				R = koo.calcStromgrenRadius(logQ, nhs_init[col])
				#R = max(R, inj[col][row] / CM2PC)
				R = koo.calcSpitzerRadius2(R, times[col]*YR2S)*CM2PC
				R = koo.calcExactIF(300.0, T, 1.0, mu, nhs_init[col], 10.0**logQ,
									times[col] * YR2S, is_spitzer=True) * CM2PC
				tableline = tableline + str(R)
			print tableline + " \\\\"

			#tableline = "& {$R_\\mathrm{stag}$}"
			#for col in range(5):
			#	tableline = tableline + " & "
			#	R = koo.calcStromgrenRadius(logQ, nhs_init[col])
			#	#R = max(R, inj[col][row] / CM2PC)
			#	R = koo.calcStagnationRadius2(R) * CM2PC
			#	tableline = tableline + str(R)
			#print tableline + " \\\\"

		if row != 8:
			print "& & & & & & \\\\"

def ratioIF_SW(isMVD):
	for row in range(9):
		mdot = mdotList[row]
		vinf = vinfList[row]

		sw_r = mvd_SW_r if isMVD else mvt_SW_r
		if_r = mvd_IF_r if isMVD else mvt_IF_r

		tableline =  str(massList[row])
		for col in range(5):
			rat = sw_r[col][row] / if_r[col][row]

			tableline = tableline + " & "
			tableline = tableline + str(rat)
		print tableline + " \\\\"

def printSlope():
	for row in range(9):
		print mvt_SW_r[:,row]
		print mvt_times
		slope, intercept, r_val, p_val, std_err = stats.linregress(np.log(mvt_SW_r[:,row]), np.log(mvt_times))
		print slope
		print std_err

#printPaperTable(mvd_IF_r, mvd_IF_t, mvd_IF_b, False, True)
#printPaperTable(mvt_IF_r, mvt_IF_t, mvt_IF_b, False, False)

#printPaperTable(mvd_SW_r, mvd_SW_t, mvd_SW_b, True, True)
#printPaperTable(mvt_SW_r, mvt_SW_t, mvt_SW_b, True, False)

#printAnalyticalSizes(mvd_SW_r, True)
#printAnalyticalSizes(mvt_SW_r, False)

#printAnalyticalTimes(True)
#printAnalyticalTimes(False)

#printCritVel(True)
#printCritVel(False)

#printAdiabaticTable(mvd_IF_r, True)
#printAdiabaticTable(mvt_IF_r, False)

#printSlope()

#ratioIF_SW(True)

windHiiVelRatio(True)
windHiiVelRatio(False)











