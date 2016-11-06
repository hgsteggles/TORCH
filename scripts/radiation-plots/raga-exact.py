import numpy as np
from scipy import integrate
from matplotlib.pylab import *

PC2CM = 3.09e18
YR2S = 3.15569e7

def spitzer_exact(t, y):
	co = 2.87e5
	c_i = 12.85e5
	T_o = 1.0e3
	T_i = 1.0e4
	mu_o = 1.0
	mu_i = 0.5
	RS = 0.314 * PC2CM

	# Output from ODE function must be a COLUMN vector, with n rows
	n = len(y)      # 1: implies its a single ODE
	dydt = np.zeros((n,1))
	dydt[0] = c_i * (RS / y[0])**(3.0 / 4.0) - c_i * (mu_i * T_o / (mu_o * T_i)) * (RS / y[0])**(- 3.0 / 4.0)
	return dydt

def raga_exact(t, y):
	co = 2.87e5
	c_i = 12.85e5
	T_o = 1.0e3
	T_i = 1.0e4
	mu_o = 1.0
	mu_i = 0.5
	RS = 0.314 * PC2CM

	inner = (4.0 / 3.0) * (RS / y)**(3.0 / 2.0) - (mu_i * T_o / (2.0 * mu_o * T_i))
	inner = max(inner, 0)

	return c_i * math.sqrt(inner)

def Euler(ta, tb, n, fa, f):
	R = [fa]
	t = [ta]
	h = (tb - ta) / float(n)
	for i in range(n):
		R.append(R[i] + h * f(t[i], R[i]))
		t.append(t[i] + h)
	return R, t

def solve_if(t_start, t_final, nsteps, R_zero, eqn):
	# Start by specifying the integrator:
	# use ``vode`` with "backward differentiation formula"
	r = integrate.ode(eqn).set_integrator('dopri5', method='bdf')

	delta_t = (t_final - t_start) / float(nsteps)

	# Number of time steps: 1 extra for initial condition
	num_steps = nsteps + 1

	# Set initial condition(s): for integrating variable and time!
	r.set_initial_value([R_zero], t_start)

	# Additional Python step: create vectors to store trajectories
	t = np.zeros((num_steps, 1))
	R = np.zeros((num_steps, 1))
	t[0] = t_start
	R[0] = R_zero

	# Integrate the ODE(s) across each delta_t timestep
	k = 1
	while r.successful() and k < num_steps:
		r.integrate(r.t + delta_t)

		# Store the results to plot later
		t[k] = r.t
		R[k] = r.y[0]
		k += 1

	return R, t

# The ``driver`` that will integrate the ODE(s):
if __name__ == '__main__':
	# Set the time range
	R_zero = 0.314 * PC2CM
	t_start = 0.0
	t_final = 3.0e6 * YR2S
	nsteps = 1000

	R_spit, t_spit = solve_if(t_start, t_final, nsteps, R_zero, spitzer_exact)
	R_raga, t_raga = Euler(t_start, t_final, nsteps, R_zero, raga_exact)

	# All done!  Plot the trajectories:
	plot_raga_i,  = plot(t_spit, R_spit, label='raga-I')
	plot_raga_ii, = plot(t_raga, R_raga, label='raga-II')

	legend(handles=[plot_raga_i, plot_raga_ii])

	grid('on')
	xlabel('Time [minutes]')
	ylabel('Concentration [mol/L]')

	show()