-- Torch Parameters
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
		setup_id       =             0,
		output_directory =           "tmp_1/",
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
		photon_rate =                2.14e43,
		wind_radius_in_cells =       10,
		mass_loss_rate =             8.42e15,
		wind_velocity =              1.13e8,
		wind_temperature =           10000,
	},
}

stars = {
{2.14e43, 8.42e15, 1.13e8},
{5.75e44, 4.41e15, 2.43e8},
{1.05e46, 5.09e16, 2.60e8},
{1.07e47, 2.73e17, 2.71e8},
{1.00e48, 1.63e18, 2.86e8},
{2.88e48, 5.15e18, 2.99e8},
{4.90e48, 9.79e18, 3.11e8},
{7.94e48, 1.88e19, 3.17e8},
{1.23e49, 3.64e19, 3.20e8},
{2.04e49, 6.79e19, 3.32e8},
{2.69e49, 1.31e20, 3.34e8},
{3.24e49, 1.84e20, 3.41e8},
{3.80e49, 2.64e20, 3.42e8},
{4.47e49, 3.69e20, 3.42e8},
{5.01e49, 4.53e20, 3.45e8},
{5.75e49, 5.63e20, 3.46e8},
{6.46e49, 7.07e20, 3.42e8}}


function param_set(n)
	local setup_id = math.floor((n-1)/17.0)  + 1
	local star_id = n - (setup_id - 1)*17

	Parameters.Integration.setup_id = setup_id
	Parameters.Integration.output_directory = "tmp_"..n.."/"
	Parameters.Star.photon_rate = stars[star_id][1]
	Parameters.Star.mass_loss_rate = stars[star_id][2]
	Parameters.Star.wind_velocity = stars[star_id][3]
end
