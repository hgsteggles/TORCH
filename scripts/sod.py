import matplotlib.pyplot as plt

def sod():
	gamma = 1.4
	rho1 = 1
	rho4 = 0.125
	u1 = 0
	u4 = 0
	p1 = 1
	p4 = 0.1

	# Calculation of flow parameters at initial condition
	a1 = (gamma*p1/rho1)**0.5
	a4 = (gamma*p4/rho4)**0.5
	M1 = u1/a1
	M4 = u4/a4

	# Secant Method for getting pressure P* 
	#p23up = (((gamma-1)/2*(u1-u4)+a1+a4)/((a1*(p1)**((2*gamma)/(gamma-1)))+(a4*(p4)**((2*gamma)/(gamma-1)))))**((2*gamma)/(gamma-1)) # upper limit
	p23up = (((gamma-1)/2*(u1-u4)+a1+a4)/((a1*(p1)**((1-gamma)/(2.0*gamma)))+(a4*(p4)**((1-gamma)/(2.0*gamma)))))**((2*gamma)/(gamma-1)) # upper limit
	p23down = (rho1*a1*p4+rho4*a4*p1-(rho1*a1*rho4*a4*(u4-u1)))/(rho1*a1+rho4*a4) # lower limit by linear theory
	s = 0
	if p23down >= p1:
		s = 1
	ss = 0
	if p23down >= p4:
		ss = 1
	if s == 1:
		m1 = rho1*a1*(1+((gamma+1)*(p23up-p1)/(2*gamma*p1)))**0.5
		m1d = rho1*a1*(1+((gamma+1)*(p23downp1)/(2*gamma*p1)))**0.5
	else:
		m1 = rho1*a1*(gamma-1)/(2*gamma)*(1-p23up/p1)/(1-(p23up/p1)**((gamma-1)/(2*gamma)))
		m1d = rho1*a1*(gamma-1)/(2*gamma)*(1-p23down/p1)/(1-(p23down/p1)**((gamma-1)/(2*gamma)))
	if ss == 1:
		m4 = rho4*a4*(1+((gamma+1)*(p23up-p4)/(2*gamma*p4)))**0.5
		m4d = rho4*a4*(1+((gamma+1)*(p23down-p4)/(2*gamma*p4)))**0.5
	else:
		m4 = rho4*a4*(gamma-1)/(2*gamma)*(1-p23up/p4)/(1-(p23up/p4)**((gamma-1)/(2*gamma)))
		m4d = rho4*a4*(gamma-1)/(2*gamma)*(1-p23down/p4)/(1-(p23down/p4)**((gamma-1)/(2*gamma)))
	p23 = (m1*p4+m4*p1-m1*m4*(u4-u1))/(m1+m4)
	f = p23up-p23
	p23d = (m1d*p4+m4d*p1-m1d*m4d*(u4-u1))/(m1d+m4d)
	ff = p23d-p23
	j = 0

	# iteration procedure starts from here
	while abs(f) > 0.000001:
		p23up = p23up-(f*(p23up-p23d)/(f-ff))
		if s == 1:
			m1 = rho1*a1*(1+((gamma+1)*(p23up-p1)/(2*gamma*p1)))**0.5
			m1d = rho1*a1*(1+((gamma+1)*(p23down-p1)/(2*gamma*p1)))**0.5
		else:
			m1 = rho1*a1*(gamma-1)/(2*gamma)*(1-p23up/p1)/(1-(p23up/p1)**((gamma-1)/(2*gamma)))
			m1d = rho1*a1*(gamma-1)/(2*gamma)*(1-p23down/p1)/(1-(p23down/p1)**((gamma-1)/(2*gamma)))
		if ss == 1:
			m4 = rho4*a4*(1+((gamma+1)*(p23up-p4)/(2*gamma*p4)))**0.5
			m4d = rho4*a4*(1+((gamma+1)*(p23down-p4)/(2*gamma*p4)))**0.5
		else:
			m4 = rho4*a4*(gamma-1)/(2*gamma)*(1-p23up/p4)/(1-(p23up/p4)**((gamma-1)/(2*gamma)))
			m4d = rho4*a4*(gamma-1)/(2*gamma)*(1-p23down/p4)/(1-(p23down/p4)**((gamma-1)/(2*gamma)))
		p23 = (m1*p4+m4*p1-m1*m4*(u4-u1))/(m1+m4)
		f = p23up-p23
		p23d = (m1d*p4+m4d*p1-m1d*m4d*(u4-u1))/(m1d+m4d)
		ff = p23d-p23
		j = j+1
		if j > 450000:
			print 'no convergance'
			break

	u23 = (m1*u1+m4*u4-(p4-p1))/(m1+m4)
	if s == 1:
		rho2 = rho1*(1+(((gamma+1)/(gamma-1))*p23/p1))/(((gamma+1)/(gamma-1))+p23/p1)
		if u23 > u1:
			print('Expansion shock-not physically possible')
			return
	else:
		rho2 = rho1*(p23/p1)**(1/gamma)
	if ss == 1:
		rho3 = rho4*(1+(((gamma+1)/(gamma-1))*p23/p4))/(((gamma+1)/(gamma-1))+p23/p4)
		if u23 < u4:
			print('Expansion shock-not physically possible')
			return
	else:
		rho3 = rho4*(p23/p4)**(1/gamma)
	a2 = (gamma*p23/rho2)**0.5
	a3 = (gamma*p23/rho3)**0.5

	# Print calculated flow quantities
	#print ('Solution of Riemann problem :-')
	#print ('P* =' + str(p23) )
	#print ('U* = ' + str(u23) )
	#print ('rho2 = ' + str(rho2) )
	#print ('rho3 = ' + str(rho3) )

	# Variable initialization
	cs12 = 0
	cs34 = 0
	expc121 = 0
	expc122 = 0
	expc341 = 0
	expc342 = 0
	# Calculation of shock/expansion speeds
	if s == 1:
		cs12l = a1*(1+(gamma+1)/(2*gamma)*(p23-p1)/p1)**0.5
		cs12 = u1-abs(cs12l)
	else:
		expc121 = u1-a1
		expc122 = u23-a2
	if ss == 1:
		cs34r = a4*(1+(gamma+1)/(2*gamma)*(p23-p4)/p4)**0.5
		cs34 = u4+abs(cs34r)
	else:
		expc341 = u4+a4
		expc342 = u23+a3

	# array construction
	maxxt = max([cs12, cs34, expc121, expc122, expc341, expc342])
	minxt = min([cs12, cs34, expc121, expc122, expc341, expc342])
	offsetxt = 0.1*(maxxt-minxt)

	xt = [0] * 6000
	u = [0] * 6000
	rho = [0] * 6000
	p = [0] * 6000
	e = [0] * 6000

	if s == 1:
		xt[0] = cs12-offsetxt
		incr = abs(offsetxt)/(1500)
		for i  in range(0, 1500):
			xt[i+1] = xt[i]+incr
			u[i] = u1
			rho[i] = rho1
			p[i] = p1
			e[i] = (p[i]/(gamma-1))/rho[i] 
		xt[1500] = cs12
		incr = abs(u23-cs12)/(1500)
		for i in range(1500, 3000):
			xt[i+1] = xt[i]+incr
			u[i] = u23
			rho[i] = rho2
			p[i] = p23
			e[i] = p[i]/(gamma-1)/rho[i]
		xt[3000] = u23
	else:
		xt[0] = expc121-offsetxt
		incr = abs(offsetxt)/(1000)
		for i in range(0, 1000):
			xt[i+1] = xt[i]+incr
			u[i] = u1
			rho[i] = rho1
			p[i] = p1
			e[i] = p[i]/(gamma-1)/rho[i]
		xt[1000] = expc121
		incr = abs(expc122-expc121)/(1000)
		for i in range(1000, 2000):
			xt[i+1] = xt[i]+incr
			if expc122 >= 0:
				u[i] = 2/(gamma+1)*(xt[i]-a1)+(gamma-1)/(gamma+1)*u1
				a = ((gamma-1)/(gamma+1)*(xt[i]-u1))+(2/(gamma+1)*a1)
			else:
				u[i] = 2/(gamma+1)*(xt[i]+a1)+(gamma-1)/(gamma+1)*u1
				a = (-(gamma-1)/(gamma+1)*(xt[i]-u1))+(2/(gamma+1)*a1)
			rho[i] = rho1*(a/a1)**(2/(gamma-1))
			p[i] = p1*(a/a1)**(2*gamma/(gamma-1))
			e[i] = p[i]/(gamma-1)/rho[i]
		xt[2000] = expc122
		incr = abs(expc122-u23)/1000
		for i in range(2000, 3000):
			xt[i+1] = xt[i]+incr
			u[i] = u23
			rho[i] = rho2
			p[i] = p23
			e[i] = p[i]/(gamma-1)/rho[i]
		xt[3000] = u23

	if ss == 1:
		incr = abs(u23-cs34)/1500
		for i in range(3000, 4500):
			xt[i] = xt[i-1]+incr
			u[i] = u23
			rho[i] = rho3
			p[i] = p23
			e[i] = p[i]/(gamma-1)/rho[i]
		xt[4500] = cs34
		incr = abs(offsetxt)/1500
		for i in range(4500, 6000):
			xt[i] = xt[i-1]+incr
			u[i] = u4
			rho[i] = rho4
			p[i] = p4
			e[i] = p[i]/(gamma-1)/rho[i]
	else:
		incr = abs(expc342-u23)/1000
		for i in range(3000, 4000):
			xt[i] = xt[i-1]+incr
			u[i] = u23
			rho[i] = rho3
			p[i] = p23
			e[i] = p[i]/(gamma-1)/rho[i]
		xt[4000] = expc342
		incr = abs(expc342-expc341)/1000
		for i in range(4000, 5000):
			xt[i] = xt[i-1]+incr
			if expc341 >= 0:
				u[i] = 2/(gamma+1)*(xt[i]-a4)+(gamma-1)/(gamma+1)*u4
				a = ((gamma-1)/(gamma+1)*(xt[i]-u4))+(2/(gamma+1)*a4)
			else:
				u[i] = 2/(gamma+1)*(xt[i]+a4)+(gamma-1)/(gamma+1)*u4
				a = (-(gamma-1)/(gamma+1)*(xt[i]-u4))+(2/(gamma+1)*a4)
			rho[i] = rho4*(a/a4)**(2/(gamma-1))
			p[i] = p4*(a/a4)**(2*gamma/(gamma-1))
			e[i] = p[i]/(gamma-1)/rho[i]
		xt[5000] = expc341
		incr = abs(offsetxt)/1000
		for i in range(5000, 6000):
			xt[i] = xt[i-1]+incr
			u[i] = u4
			rho[i] = rho4
			p[i] = p4
			e[i] = p[i]/(gamma-1)/rho[i]

	x = [0] * 6000
	for i  in range(0, 6000):
		x[i] = 0.5 + 0.25*xt[i];
	x[0] = 0
	return (x, rho)

#def main():
#		figformat = 'png'
#		outputfile='sod.png'
#		sod_x = [0] * 6000
#		sod_rho = [0] * 6000
#		sod_x, sod_rho = sod()
#		for i in range(0, 6000):
#			sod_x[i] = 100*sod_x[i]
#		plt.plot(sod_x, sod_rho)
#		plt.savefig(outputfile, format=figformat)
#
#main()

