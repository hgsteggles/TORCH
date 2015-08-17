hydrogenMass = 1.674e-24
specificGasConstant = 8.314462e7
PC2CM = 3.09e18
H = 0.05*PC2CM
T = 300
nHI = 8000

function initialise(x, y, z, xs, ys, zs)
	local z = (y - ys)/H
	local den = nHI*hydrogenMass*math.exp(z)
	local pre = specificGasConstant*den*T
	local hii = 0
	local v0 = 0
	local v1 = 0
	local v2 = 0
	return den, pre, hii, v0, v1, v2
end
