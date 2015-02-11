hydrogenMass = 1.674e-24
specificGasConstant = 8.314462e7
G = 6.6726e-8
PC2CM = 3.09e18
H = 0.05*PC2CM
T = 300
nHI = 80000

alpha = 1
rc = 0.01*PC2CM
rc2 = rc*rc
RS = 0.35*PC2CM

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
	return den, pre, hii, v0, v1, v2
end

nHI_array = {10000, 45000, 80000}
RS_array = {}
for i=1, 10 do
	RS_array[i] = (0.35 + 0.01*(i-1))*PC2CM
end

function setup_set(n)
	local nHI_id = math.floor((n-1)/10.0)  + 1
	local RS_id = n - (nHI_id - 1)*10

	nHI = nHI_array[nHI_id]
	RS = RS_array[RS_id]

	n0 = nHI*math.pow(1 + RS*RS/rc2, alpha)
	pre0 = specificGasConstant*nHI*hydrogenMass*T
end
