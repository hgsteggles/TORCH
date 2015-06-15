import numpy as np

def write_file(filename, nHI, alpha, RS):
	f = open(filename, 'w')
	
	string ="""hydrogenMass = 1.674e-24
specificGasConstant = 8.314462e7
G = 6.6726e-8
PC2CM = 3.09e18
H = 0.05*PC2CM
T = 300
nHI = """ + str(nHI) + "\n" + """
alpha = """ + str(alpha) + "\n" + """
rc = 0.01*PC2CM
rc2 = rc*rc
RS = """ + str(RS) + "*PC2CM\n" + """		
n0 = nHI*math.pow(1 + RS*RS/rc2, alpha)

pre0 = specificGasConstant*nHI*hydrogenMass*T

function initialise(x, y, z, xs, ys, zs)
	local dy = RS + (ys - y)
	local R2 = x*x + dy*dy
	local R = math.sqrt(R2)

	local den = n0*hydrogenMass*math.pow(1 + R2/rc2, -alpha)
	local pre = pre0
	local hii = 0
	local v0 = 0
	local v1 = 0
	local v2 = 0
	local grav0 = 0
	local grav1 = 0
	local grav2 = 0
	return den, pre, hii, v0, v1, v2, grav0, grav1, grav2
end """
	f.write(string)

nHI_array = [10000, 40000, 80000]
RS_array = np.arange(0.35, 0.44, 0.01)
alpha = 2

i = 1

for nHI in nHI_array:
	for RS in RS_array:
		write_file("setup_"+str(i)+".lua", nHI, alpha, RS)
		i += 1;


