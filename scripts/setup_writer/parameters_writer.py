import numpy as np
import math

def write_file(i, setup_i, photon_rate, mass_loss_rate, wind_velocity):
	f = open("parameters_"+str(i)+".lua", 'w')
	string = """-- Torch Parameters
Parameters = {
	Integration = {
		ambient_no_density =         8.0e3,
		density_scale =              1.0e-20,
		pressure_scale =             1.0e-28,
		time_scale =                 1.0e11,
		spatial_order =              1,
		temporal_order =             2,
		simulation_time =            2.0*3.15569e12,
		radiation_on =               true,
		cooling_on =                 true,
		debug =                      true,
		lua_setup_file =             "setup_""" + str(setup_i) + """.lua",
		output_directory =           "tmp_""" + str(i) + """/",
		initial_conditions =         "",
	},
	Grid = {
		no_dimensions =              2,
		no_cells_x =                 150,
		no_cells_y =                 200,
		no_cells_z =                 1,
		geometry =                   "cylindrical",
		side_length =                1.854e18,
		left_boundary_condition_x =  "reflecting",
		left_boundary_condition_y =  "outflow",
		left_boundary_condition_z =  "reflecting",
		right_boundary_condition_x = "outflow",
		right_boundary_condition_y = "reflecting",
		right_boundary_condition_z = "outflow",
	},
	Hydrodynamics = {
		gamma =                      1.67,
		density_floor =              1.0e-30,
		pressure_floor =             1.0e-22,
		temperature_floor =          0.1,
		riemann_solver =             "RotatedHLLC",
		slope_limiter =              "albada",
	},
	Radiation = {
		K1 =                         0.2,
		K2 =                         0.0,
		K3 =                         0.0,
		K4 =                         0.0,
		photoion_cross_section =     6.3e-18,
		case_b_recombination_coeff = 2.59e-13,
		tau_0 =                      0.6,
		minimum_hii_fraction =       0,
		temperature_hi =             300,
		temperature_hii =            10000,
		mass_fraction_hydrogen =     1.0,
		integration_scheme =         "implicit",
		collisions_on =              false,
		coupling =                   "neq",
	},
	Thermodynamics = {
		thermo_hii_switch =          1.0e-2,
		heating_amplification =      1,
	},
	Star = {
		on =                         true,
		cell_position_x =            0,
		cell_position_y =            120,
		cell_position_z =            0,
		snap_to_face_left_x =        true,
		snap_to_face_left_y =        false,
		snap_to_face_left_z =        false,
		snap_to_face_right_x =       false,
		snap_to_face_right_y =       false,
		snap_to_face_right_z =       false,
		photon_energy =              2.976e-11,
		photon_rate =                """+str(photon_rate)+""",
		wind_radius_in_cells =       10,
		mass_loss_rate =             """+str(mass_loss_rate)+""",
		wind_velocity =              """+str('%.2e' % wind_velocity)+""",
		wind_temperature =           10000,
	},
}"""
	f.write(string)

round_to_n = lambda x, n: round(x, -int(math.floor(math.log10(x))) + (n - 1))

data = np.array([[43.33, 8.42e15, 1.13e8],
[44.76, 4.41e15, 2.43e8],
[46.02, 5.09e16, 2.60e8],
[47.03, 2.73e17, 2.71e8],
[48.00, 1.63e18, 2.86e8],
[48.46, 5.15e18, 2.99e8],
[48.69, 9.79e18, 3.11e8],
[48.90, 1.88e19, 3.17e8],
[49.09, 3.64e19, 3.20e8],
[49.31, 6.79e19, 3.32e8],
[49.43, 1.31e20, 3.34e8],
[49.51, 1.84e20, 3.41e8],
[49.58, 2.64e20, 3.42e8],
[49.65, 3.69e20, 3.42e8],
[49.70, 4.53e20, 3.45e8],
[49.76, 5.63e20, 3.46e8],
[49.81, 7.07e20, 3.42e8]])

for datum in data:
	print round_to_n(math.pow(10, datum[0]), 3)

