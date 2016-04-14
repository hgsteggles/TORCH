import numpy as np
import math

print_sizes = False
print_mvd_table = False
print_mvt_table = True

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

nhs = [0.8e4, 1.6e4, 3.2e4, 6.4e4, 12.8e4]
times = [2.0e4, 4.0e4, 6.0e4, 8.0e4, 10.0e4]

def calcBubbleRadius(mdot, vinf, rho_0, time):
	return 0.88*pow(mdot * vinf * vinf * time * time * time / rho_0, 0.2)

def calcBubbleRadiusRadiative(mdot, vinf, rho_0, time):
	return pow(3.0 * 0.5 * mdot * vinf * time * time / (math.pi * rho_0), 0.25)

def calcBubbleContact(mdot, vinf, rho_0, time):
	rho = 0.1*rho_0
	T = 8000.0
	return pow(3.0 * 0.5 * mdot * vinf * vinf * time / (10.0 * math.pi * (R_GAS * T * rho / 0.5)), 1.0/3.0)

def calcIFrontRadius(logQ, nh):
	T = 8000
	alphaB = 2.59e-13*math.pow(T/10000.0, -0.7)
	return pow(3.0 * pow(10.0, logQ) / (4.0 * math.pi * nh * nh * alphaB), 1.0/3.0)

def calcSpitzerRadius(logQ, nh, t):
	RS = calcIFrontRadius(logQ, nh)
	cs = math.sqrt(2.0*R_GAS*8000.0)
	ts = RS/cs
	return RS*math.pow(1.0 + ((7.0*t)/(4.0*ts)), 4.0/7.0)

def calcCorrectedSpitzer(mdot, t, nh, R):
	V = (4.0 / 3.0) * math.pi * R * R * R
	return math.pow(1.0 + (mdot * t/(mh * nh * V)), -2.0/3.0)*R

def calcPressureEqual(mdot, vinf, rho_0):
	T = 8000.0
	#rho = 0.1 * rho_0
	rho = rho_0
	return math.sqrt(3.0 * 0.5 * mdot * vinf * 0.5 * 0.5 / (16.0 * math.pi * rho * R_GAS * R_GAS * T * T))

def calcCritVel(mdot, vinf, rho_0):
	#rho = 0.1 * rho_0
	rho = rho_0
	C1 = 6.0e-35
	return math.pow(rho * 0.5 * mdot * vinf * vinf / (2.0 * math.pi * C1 * C1), 1.0/11.0)

def calcPressureEqual2(mdot, vinf, rho_0):
	T = 8000.0
	#rho = 0.1 * rho_0
	rho = rho_0
	c5 = math.pow(math.sqrt(R_GAS * T / 0.5), 5.0)
	return 0.142 * math.sqrt(0.5 * mdot * vinf * vinf / (rho * c5))

def finalRadius(mdot, vinf, rho_0):
	tp25 = math.pow(calcPressureEqual2(mdot, vinf, rho_0), 2.0/5.0)
	return 0.92 * math.pow(0.5 * mdot * vinf * vinf / (rho_0 * math.pow(vinf, 5.0/3.0)), 0.3) * tp25

######### MVD

#### High Temp Wind

mvd_xrange = []
mvd_xrange.append([0.10, 0.24, 0.46, 0.70, 0.90, 1.0, 1, 1, 1])
mvd_xrange.append([0.06, 0.15, 0.38, 0.60, 0.80, 0.9, 1, 1, 1])
mvd_xrange.append([0.04, 0.11, 0.25, 0.45, 0.65, 0.9, 1, 1, 1])
mvd_xrange.append([0.03, 0.06, 0.20, 0.35, 0.55, 0.8, 1, 1, 1])
mvd_xrange.append([0.02, 0.05, 0.14, 0.25, 0.45, 0.6, 0.8, 1, 1])
mvd_xrange = np.array(mvd_xrange, np.float64)

mvd_ref = 33.6

mvd_hsize = []
mvd_hsize.append([16.5, 4.25, 5.50, 5.00, 6.00, None, None, None, None])
mvd_hsize.append([19.5, 8.25, 4.00, 5.00, 6.25, 10.5, None, None, None])
mvd_hsize.append([16.5, 7.00, 7.50, 6.00, 7.00, 11.0, 18.5, None, None])
mvd_hsize.append([15.0, 10.5, 6.75, 6.00, 6.50, 9.50, 12.2, None, None])
mvd_hsize.append([14.5, 5.00, 8.50, 5.50, 7.50, 16.0, 15.0, 25.0, None])
mvd_hsize = np.array(mvd_hsize, np.float64)

mvd_vup = []
mvd_vup.append([6.00, 2.25, 3.00, 1.50, 1.50, None, None, None, None])
mvd_vup.append([5.50, 5.50, 1.00, 2.50, 2.00, 3.00, None, None, None])
mvd_vup.append([5.00, 3.00, 3.00, 3.00, 3.00, 3.00, 5.00, None, None])
mvd_vup.append([7.00, 2.00, 3.50, 2.50, 4.00, 3.00, 4.00, None, None])
mvd_vup.append([6.00, 2.00, 3.00, 3.00, 3.25, 4.00, 4.00, 6.00, None])
mvd_vup = np.array(mvd_vup, np.float64)

mvd_vdown = []
mvd_vdown.append([11.0, 12.0, 10.5, 10.5, 13.0, None, None, None, None])
mvd_vdown.append([10.0, 4.00, 10.0, 9.00, 11.0, 14.5, None, None, None])
mvd_vdown.append([10.0, 7.00, 8.50, 9.00, 10.0, 13.5, None, None, None])
mvd_vdown.append([10.5, 8.00, 6.50, 8.00, 8.50, 10.5, 13.5, None, None])
mvd_vdown.append([10.5, 4.50, 5.00, 9.00, 8.50, 9.50, 12.5, None, None])
mvd_vdown = np.array(mvd_vdown, np.float64)


### Ionisation Front

mvd_i_hsize = []
mvd_i_hsize.append([25.0, 24.0, 25.0, 25.0, 30.0, None, None, None, None])
mvd_i_hsize.append([26.0, 27.0, 22.0, 23.0, 27.0, 32.5, None, None, None])
mvd_i_hsize.append([23.0, 24.0, 24.0, 23.0, 25.5, 26.0, 30.0, None, None])
mvd_i_hsize.append([23.5, 26.5, 21.0, 21.5, 23.5, 22.5, 23.0, None, None])
mvd_i_hsize.append([21.0, 19.5, 19.5, 21.5, 22.0, 23.5, 22.5, 26.5, None])
mvd_i_hsize = np.array(mvd_i_hsize, np.float64)

mvd_i_vup = []
mvd_i_vup.append([11.5, 11.0, 11.0, 9.50, 10.0, None, None, None, None])
mvd_i_vup.append([10.5, 12.0, 8.50, 9.00, 10.0, 10.0, None, None, None])
mvd_i_vup.append([8.50, 12.0, 10.5, 10.5, 10.0, 8.5, 9.00, None, None])
mvd_i_vup.append([11.0, 12.0, 10.5, 9.00, 10.5, 9.00, 8.00, None, None])
mvd_i_vup.append([9.50, 9.50, 8.50, 9.50, 10.0, 10.0, 8.50, 8.00, None])
mvd_i_vup = np.array(mvd_i_vup, np.float64)

mvd_i_vdown = []
mvd_i_vdown.append([14.5, 15.0, 14.0, 15.0, None, None, None, None, None])
mvd_i_vdown.append([14.5, 16.0, 13.5, 13.5, None, None, None, None, None])
mvd_i_vdown.append([13.5, 13.0, 15.0, 13.5, 15.5, None, None, None, None])
mvd_i_vdown.append([12.0, 14.0, 12.0, 13.0, 13.5, 14.5, None, None, None])
mvd_i_vdown.append([13.0, 10.5, 11.0, 13.0, 12.5, 13.5, 14.0, None, None])
mvd_i_vdown = np.array(mvd_i_vdown, np.float64)


######### MVT

### High Temp Wind

mvt_xrange = []
mvt_xrange.append([0.040, 0.10, 0.20, 0.28, 0.40, 0.52, 0.7, 1.0, 1.0])
mvt_xrange.append([0.035, 0.09, 0.22, 0.32, 0.50, 0.80, 1.0, 1.0, 1.0])
mvt_xrange.append([0.040, 0.10, 0.24, 0.40, 0.65, 0.90, 1.0, 1.0, 1.0])
mvt_xrange.append([0.040, 0.08, 0.25, 0.48, 0.74, 1.00, 1.0, 1.0, 1.0])
mvt_xrange.append([0.040, 0.10, 0.25, 0.54, 0.80, 1.00, 1.0, 1.0, 1.0])
mvt_xrange = np.array(mvt_xrange, np.float64)


mvt_ref = 34.0

mvt_hsize = []
mvt_hsize.append([19.0, 7.00, 5.50, 8.00, 7.00, 11.5, 20.0, 22.0, None])
mvt_hsize.append([20.0, 6.50, 5.50, 8.00, 8.00, 8.50, 12.0, None, None])
mvt_hsize.append([17.5, 8.00, 6.00, 7.00, 8.50, 9.00, 14.0, None, None])
mvt_hsize.append([19.0, 7.00, 7.50, 5.50, 6.00, 8.00, None, None, None])
mvt_hsize.append([19.0, 6.00, 7.50, 5.50, 5.00, None, None, None, None])
mvt_hsize = np.array(mvt_hsize, np.float64)

mvt_vup = []
mvt_vup.append([7.50, 3.00, 3.00, 5.00, 3.00, 6.00, 6.00, 6.00, None])
mvt_vup.append([6.00, 4.00, 2.00, 4.00, 4.00, 3.00, 5.00, None, None])
mvt_vup.append([5.00, 5.00, 4.00, 3.00, 3.00, 3.00, 5.00, None, None])
mvt_vup.append([8.00, 4.50, 3.00, 3.00, 2.00, 3.00, None, None, None])
mvt_vup.append([8.00, 6.00, 6.00, 3.00, 2.00, None, None, None, None])
mvt_vup = np.array(mvt_vup, np.float64)

mvt_vdown = []
mvt_vdown.append([11.0, 7.50, 8.50, 6.00, 9.00, 12.0, 13.5, None, None])
mvt_vdown.append([9.50, 10.0, 12.0, 10.5, 12.5, 14.0, None, None, None])
mvt_vdown.append([11.0, 7.50, 7.50, 12.0, 12.0, 14.0, None, None, None])
mvt_vdown.append([10.5, 3.00, 11.0, 12.5, 13.0, None, None, None, None])
mvt_vdown.append([11.0, 2.00, 9.00, 12.0, 14.0, None, None, None, None])
mvt_vdown = np.array(mvt_vdown, np.float64)


### Ionisation Front

mvt_i_hsize = []
mvt_i_hsize.append([27.0, 24.5, 23.5, 27.0, 26.0, 27.0, 25.0, 24.5, None])
mvt_i_hsize.append([28.0, 31.0, 27.0, 31.5, 30.0, 25.5, 16.5, None, None])
mvt_i_hsize.append([25.5, 27.5, 27.5, 29.0, 27.5, 27.5, 32.0, None, None])
mvt_i_hsize.append([28.0, 33.0, 28.0, 27.0, 28.0, 29.0, None, None, None])
mvt_i_hsize.append([27.5, 26.0, 29.0, 26.0, 28.0, None, None, None, None])
mvt_i_hsize = np.array(mvt_i_hsize, np.float64)

mvt_i_vup = []
mvt_i_vup.append([12.5, 13.0, 10.0, 13.0, 11.0, 11.0, 9.00, 8.00, None])
mvt_i_vup.append([12.0, 15.5, 12.0, 14.0, 12.0, 9.00, 8.00, None, None])
mvt_i_vup.append([10.0, 14.0, 12.0, 13.0, 11.0, 9.00, 9.00, None, None])
mvt_i_vup.append([13.5, 16.0, 12.5, 12.0, 11.0, 9.00, None, None, None])
mvt_i_vup.append([12.0, 13.0, 13.0, 11.5, 11.0, None, None, None, None])
mvt_i_vup = np.array(mvt_i_vup, np.float64)

mvt_i_vdown = []
mvt_i_vdown.append([16.0, 14.0, 13.5, 14.0, 14.5, 16.0, 15.0, None, None])
mvt_i_vdown.append([14.5, 16.5, 16.5, 17.0, 16.5, 16.0, None, None, None])
mvt_i_vdown.append([16.5, 14.5, 16.5, 16.5, 16.0, None, None, None, None])
mvt_i_vdown.append([16.5, 16.5, 16.5, 16.0, 16.0, None, None, None, None])
mvt_i_vdown.append([15.5, 17.0, 17.0, 15.5, None, None, None, None, None])
mvt_i_vdown = np.array(mvt_i_vdown, np.float64)

mvt_param_index = range(19, 28, 1)
mvt_snapshots = [10, 20, 30, 40, 50]

print type(mvd_xrange[0][0])

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

if print_sizes:
	for index in range(45):
		col = int(index / 9)
		row = index % 9

		SW_l = 'None' if  mvd_hsize[col][row] == None \
			else '{:.5f}'.format(xrange[col][row] *
								 0.5 * mvd_hsize[col][row] / mvd_ref)
		SW_r = SW_l
		SW_b = 'None' if  mvd_vdown[col][row] == None \
			else '{:.5f}'.format(xrange[col][row] * mvd_vdown[col][row] / mvd_ref)
		SW_t = 'None' if  mvd_vup[col][row] == None \
			else '{:.5f}'.format(xrange[col][row] * mvd_vup[col][row] / mvd_ref)

		IF_l = 'None' if  mvd_i_hsize[col][row] == None \
			else '{:.5f}'.format(xrange[col][row] *
								 0.5 * mvd_i_hsize[col][row] / mvd_ref)
		IF_r = IF_l
		IF_b = 'None' if  mvd_i_vdown[col][row] == None \
			else '{:.5f}'.format(xrange[col][row] *
								 mvd_i_vdown[col][row] / mvd_ref)
		IF_t = 'None' if  mvd_i_vup[col][row] == None \
			else '{:.5f}'.format(xrange[col][row] *
								 mvd_i_vup[col][row] / mvd_ref)

		print str(index + 1) + ",25," +  \
			  SW_l + "," + SW_r + "," + SW_b + "," + SW_t + "," + \
			  IF_l + "," + IF_r + "," + IF_b + "," + IF_t

	for col in range(5):
		for row in range(9):
			SW_l = 'None' if  mvt_hsize[col][row] == None \
				else '{:.5f}'.format(xrange[col][row] *
									 0.5 * mvt_hsize[col][row] / mvt_ref)
			SW_r = SW_l
			SW_b = 'None' if  mvt_vdown[col][row] == None \
				else '{:.5f}'.format(xrange[col][row] *
									 mvt_vdown[col][row] / mvt_ref)
			SW_t = 'None' if  mvt_vup[col][row] == None \
				else '{:.5f}'.format(xrange[col][row] *
									 mvt_vup[col][row] / mvt_ref)

			IF_l = 'None' if  mvt_i_hsize[col][row] == None \
				else '{:.5f}'.format(xrange[col][row] *
									 0.5 * mvt_i_hsize[col][row] / mvt_ref)
			IF_r = IF_l
			IF_b = 'None' if  mvt_i_vdown[col][row] == None \
				else '{:.5f}'.format(xrange[col][row] *
									 mvt_i_vdown[col][row] / mvt_ref)
			IF_t = 'None' if  mvt_i_vup[col][row] == None \
				else '{:.5f}'.format(xrange[col][row] *
									 mvt_i_vup[col][row] / mvt_ref)

			print str(mvt_param_index[row]) + "," + str(mvt_snapshots[col]) + "," +  \
				  SW_l + "," + SW_r + "," + SW_b + "," + SW_t + "," + \
				  IF_l + "," + IF_r + "," + IF_b + "," + IF_t

def printTable(r, t, b, isWind, isMVD):
	for row in range(9):
		tableline = "& {$R_\\mathrm{\\theta = 0}\ /\ \si{\parsec}$}"
		for col in range(5):
			tableline = tableline + " & "
			if t[col][row] != None and not np.isnan(t[col][row]):
				tableline = tableline + str(t[col][row])
			else:
				tableline = tableline + "{\\textemdash}"

		print tableline + " \\\\"

		tableline =  "\multirow{4}{*}{\\num{" + str(massList[row]) + "}} & {$R_\\mathrm{\\theta = \\frac{\\pi}{2}}\ /\ \si{\parsec}$}"
		for col in range(5):
			tableline = tableline + " & "
			if r[col][row] != None and not np.isnan(r[col][row]):
				tableline = tableline + str(r[col][row])
			else:
				tableline = tableline + "{\\textemdash}"

		print tableline + " \\\\"

		tableline = "& {$R_\\mathrm{\\theta = \\pi}\ /\ \si{\parsec}$}"
		for col in range(5):
			tableline = tableline + " & "
			if b[col][row] != None and not np.isnan(b[col][row]):
				tableline = tableline + str(b[col][row])
			else:
				tableline = tableline + "{\\textemdash}"

		print tableline + " \\\\"

		if isWind:
			#tableline = "& {$R_\\mathrm{shock}\ /\ \si{\parsec}$}"
			#for col in range(5):
			#	tableline = tableline + " & "
			#	R = 0
			#	if isMVD:
			#		R = calcBubbleContact(mdotList[row], vinfList[row], mh*nhs[col], 5.0e4*YR2S)*CM2PC
			#		#R = calcBubbleRadius(mdotList[row], vinfList[row], mh*nhs[col], 5.0e4*YR2S)*CM2PC
			#		#R = calcBubbleRadiusRadiative(mdotList[row], vinfList[row], mh*nhs[col], 5.0e4*YR2S)*CM2PC
			#	else:
			#		R = calcBubbleContact(mdotList[row], vinfList[row], mh*nhs[2], times[col]*YR2S)*CM2PC
			#		#R = calcBubbleRadius(mdotList[row], vinfList[row], mh*nhs[2], times[col]*YR2S)*CM2PC
			#		#R = calcBubbleRadiusRadiative(mdotList[row], vinfList[row], mh*nhs[2], times[col]*YR2S)*CM2PC
			#
			#	tableline = tableline + str(R)

			#print tableline + " \\\\"

			#tableline = "& {$v_\\mathrm{cr}\ /\ \si{\\kilo\\meter\\per\\second}$}"
			#for col in range(5):
			#	tableline = tableline + " & "
			#	R = 0
			#	if isMVD:
			#		R = calcCritVel(mdotList[row], vinfList[row], mh*nhs[col])/100000.0
			#	else:
			#		R = calcCritVel(mdotList[row], vinfList[row], mh*nhs[2])/100000.0
			#	tableline = tableline + str(R)

			#print tableline + " \\\\"

			tableline = "& {$R_\\mathrm{P}\ /\ \\si{\\parsec}$}"
			for col in range(5):
				tableline = tableline + " & "
				R = 0
				if isMVD:
					R = finalRadius(mdotList[row], vinfList[row], mh*nhs[col]) * CM2PC
				else:
					R = finalRadius(mdotList[row], vinfList[row], mh*nhs[2]) * CM2PC
				tableline = tableline + str(R)

			print tableline + " \\\\"

		if not isWind:
			tableline = "& {$R_\\mathrm{i}\ /\ \si{\parsec}$}"
			for col in range(5):
				tableline = tableline + " & "
				R = 0
				if isMVD:
					R = calcSpitzerRadius(logQList[row], nhs[col], 5.0e4*YR2S)*CM2PC
				else:
					R = calcSpitzerRadius(logQList[row], nhs[2], times[col]*YR2S)*CM2PC
				tableline = tableline + str(R)

			print tableline + " \\\\"

			tableline = "& {$R_\\mathrm{S}\ /\ \si{\parsec}$}"
			for col in range(5):
				tableline = tableline + " & "
				R = 0
				if isMVD:
					R = calcIFrontRadius(logQList[row], nhs[col])*CM2PC
				else:
					R = calcIFrontRadius(logQList[row], nhs[2])*CM2PC
				tableline = tableline + str(R)

			print tableline + " \\\\"

		if row != 8:
			print "& & & & & & \\\\"

#printTable(mvd_SW_r, mvd_SW_t, mvd_SW_b, True, True)
#printTable(mvt_SW_r, mvt_SW_t, mvt_SW_b, True, False)
#printTable(mvd_IF_r, mvd_IF_t, mvd_IF_b, False, True)
printTable(mvt_IF_r, mvt_IF_t, mvt_IF_b, False, False)


















