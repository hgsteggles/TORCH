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
	local pre = specificGasConstant*den*T

	local hii = 0
	local v0 = 0
	local v1 = 0
	local v2 = 0


	local g_coeff = -2.0*n0*hydrogenMass*specificGasConstant*T*alpha/rc2
	local g_r = g_coeff*R*math.pow(1 + R2/rc2, -alpha-1.0)

	local grav0 = g_r*x/R
	local grav1 = -g_r*dy/R
	local grav2 = 0

	return den, pre, hii, v0, v1, v2, grav0, grav1, grav2
end
