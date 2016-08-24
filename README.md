# TORCH
### The Operator-split Radiation, Cooling, Hydro-code

******************************

####Overview
TORCH is a 3D Eulerian fixed grid fluid dynamics code. The grid is a collection of finite  elements, called grid cells, that each hold fluid state information. The hydrodynamics are solved using a rotated hybrid HLLC-HLL Riemann solver ([Nishikawa & Kitamura 2008](#N8)) to calculate fluxes on each grid cell face. Ionisation from point source radiation is implicitly solved and the column densities required for  this are calculated via an interpolative ray tracing scheme ([Mellema et al. 2006](#M6)). Heating/cooling from atomic processes is calculated using the approximate functions in [Henney et al. (2009)](#H9).

####Example Usage
First of all, the directory tree should look like this:
```
.
├── torch
├── refdata
|   ├── parameters.lua
|   └── setup.lua
```

To run, execute torch in an mpi environment (or not). For example, using 8 logical cores:
```bash
mpirun -np 8 ./torch
```
By default TORCH reads in the configuration files, "refdata/parameters.lua" and "refdata/setup.lua".
You can specify your own configuration files:
```bash
mpirun -np 8 ./torch --paramfile=/path/to/parameters.lua --setupfile=/path/to/setup.lua
```
In the configuration files, cgs units are assumed.  

#####Setup

For example, to set up a 2D cylindrically symmetric 150x200 mesh with a star located at grid coordinates (0, 110) parameters could be:

```lua
-- refdata/parameters.lua

PC2CM = 3.09e18
YR2S  = 3.15569e7

Parameters = {
	Integration = {
		density_scale =              1.0e-20,
		pressure_scale =             1.0e-28,
		time_scale =                 1.0e11,
		spatial_order =              1,
		temporal_order =             2,
		simulation_time =            2e5 * YR2S,
		radiation_on =               true,
		cooling_on =                 true,
		debug =                      false,
		output_directory =           "tmp",
		initial_conditions =         "",
	},
	Grid = {
		no_dimensions =              2,
		no_cells_x =                 150,
		no_cells_y =                 200,
		no_cells_z =                 1,
		geometry =                   "cylindrical",
		side_length =                0.5 * PC2CM,
		left_boundary_condition_x =  "reflecting",
		left_boundary_condition_y =  "outflow",
		left_boundary_condition_z =  "reflecting",
		right_boundary_condition_x = "outflow",
		right_boundary_condition_y = "outflow",
		right_boundary_condition_z = "outlfow",
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
		heating_amplification =      1,
		thermo_hii_switch =          1e-2,
		thermo_subcycling =          true,
		min_temp_initial_state =     true,
	},
	Star = {
		on =                         true,
		cell_position_x =            0,
		cell_position_y =            110,
		cell_position_z =            0,
		snap_to_face_left_x =        true,
		snap_to_face_left_y =        false,
		snap_to_face_left_z =        false,
		snap_to_face_right_x =       false,
		snap_to_face_right_y =       false,
		snap_to_face_right_z =       false,
		photon_energy =              2.976e-11,
		photon_rate =                4.9e+48,
		wind_radius_in_cells =       10,
		mass_loss_rate =             9.79e+18,
		wind_velocity =              311000000.0,
		wind_temperature =           10000,
	},
}
```
The simulated span of time in this example is 20,000 years. Snapshots of the star  would be taken at equally spaced intervals during this time and stored in the tmp  directory.
 
An example of how to set up a problem is given in "refdata/setup.lua". For a star  offset in a spherically symmetric density gradient:

```lua
-- refdata/setup.lua

hydrogenMass = 1.674e-24
specificGasConstant = 8.314462e7
PC2CM = 3.09e18

H = 0.05 * PC2CM
T = 300.0
nHI = 32000.0

alpha = 1
rc = 0.01 * PC2CM
rc2 = rc * rc
RS = 0.35 * PC2CM

n0 = nHI * math.pow(1 + RS * RS / rc2, alpha)

pre0 = specificGasConstant * nHI * hydrogenMass * T

function initialise(x, y, z, xs, ys, zs)
	local dy = RS + (ys - y)
	local R2 = x * x + dy * dy
	local R = math.sqrt(R2)

	local den = n0 * hydrogenMass * math.pow(1 + R2 / rc2, -alpha)
	local pre = pre0

	local hii = 0
	local v0 = 0
	local v1 = 0
	local v2 = 0

	local grav0 = 0
	local grav1 = 0
	local grav2 = 0

	return den, pre, hii, v0, v1, v2, grav0, grav1, grav2
end
```

The function ```initialise``` takes in six arguments: x (or r-polar), y (or z-polar) and z coordinates of the cell and the coordinates of the star (all in cm). Density, pressure, ionised hydrogen fraction, velocity components, and gravitational acceleration components are returned (in cgs units). This script is executed by Torch in order to set up the fluid variables in a grid.

#####Output
Torch outputs compressed data files in a specified directory (```output_directory```). The header contains 4 lines; the first line is the simulation time in seconds and the next three lines give the number of grid cells along the x, y and z directions of the mesh. After the header, grid cell data is displayed in columns. The first ND columns are the position coordinates of the grid cell, where ND is the number of dimensions. Next is density, pressure and HII fraction. Then the last ND columns are the fluid velocity components. All output is in cgs units.

After 5,000 years the solution to the setup given above looks like this:

![SolutionImage](four-panel-d24-t025.png)
\[Image produced using matplotlib.\]

####Getting Started

TODO

####Advanced Usage

The parameters not included in this table should not be modified. Asterisks are wildcard strings.

#####Basic
| Parameter                     | Notes                                     |
| :---------------------------- | :---------------------------------------- |
| ```*_scale```                 | Chosen such that code units of order of unity. |
| ```radiation_on```            | Simulate radiative transfer? |
| ```cooling_on```              | Simulate heating and cooling? |
| ```simulation_time```         | Span of time in seconds over which you want to simulate the fluid. |
| ```output_directory```        | Directory for output data. |
| ```initial_conditions```      | Data file to read a problem setup. Set to empty string to use setup.lua config.|
| ```no_dimensions```           | No. of dimensions in numerical grid. |
| ```no_cells_x```              | No. of cells along the x (or polar r) axis. |
| ```no_cells_y```              | No. of cells along the y (or polar z) axis. |
| ```no_cells_z```              | No. of cells along the z axis. |
| ```geometry```                | cartesian, cylindrical or spherical. |
| ```side_length```             | Length, in cm, of the x (or r) axis. Note: cells are cubic.|
| ```*_boundary_condition_*```  | reflecting, free, inflow, outflow or periodic. |
| ```gamma```                   | Heat capacity ratio. |
| ```*_floor```                 | Minimum values in cgs units. Must be positive and non-zero. |
| ```mass_fraction_hydrogen```  | Fraction, by mass, of gas in a cell that is hydrogen. |
| ```collisions_on```           | Include collisional ionisations. |
| ```heating_amplification```   | Fraction of calculated heating/cooling that is injected into fluid. |
| ```on```                      | Include star? |
| ```cell_position_x```         | Star's position along x-axis in grid coordinates. |
| ```cell_position_y```         | Star's position along y-axis in grid coordinates. |
| ```cell_position_z```         | Star's position along z-axis in grid coordinates. |
| ```snap_to_face_*```          | Normally star sits in centre of specified cell. Snap to left if that boundary is reflecting. Snap to left in x direction if cylindrical or spherical geometries. |
| ```photon_energy```           | Energy of each photon emitted by star. |
| ```photon_rate```             | Rate of photons emitted by star. |
| ```wind_radius_in_cells```    | Radius within which to inject stellar wind energy. Should be > 10 cells in 2 or 3 dimensions so that wind region is roughly spherical. |
| ```mass_loss_rate```          | Mass loss rate of star. |
| ```wind_velocity```           | Terminal velocity of the stellar wind. |
| ```wind_temperature```        | Temperature of the stellar wind region. |

##### Advanced
| Parameter                     | Notes                                     |
| :---------------------------- | :---------------------------------------- |
| ```spatial_order```           | The order of spatial reconstruction. No reconstruction with 0 and linear reconstruction with 1. |
| ```debug_on```                | Output debugging info to console? |
| ```riemann_solver```          | HLL, HLLC or RotatedHLLC. |
| ```slope_limiter```           | albada, superbee, monotonised_central, minmod or maxmod. |
| ```integration_scheme```      | Radiation integration scheme. implicit or explicit. |
| ```coupling```                | Coupling between radiation and hydrodynamics. neq (Non-equilibrium) or tti (two-temperature isothermal). |

####Goals
* Output in HDF5 data format.
* Include molecular hydrogen.

####Developer info

Harrison Steggles, University of Leeds (PhD student).

####References
<a name="H9"></a>Henney, W. J., Arthur, S. J., de Colle, F., & Mellema, G. 2009, MNRAS, 398, 157 ([link](http://mnras.oxfordjournals.org/content/398/1/157.full.pdf+html))  
<a name="M6"></a>Mellema, G., Iliev, I. T., Alvarez, M. A., & Shapiro, P. R. 2006, New A, 11, 374 ([link](http://arxiv.org/pdf/astro-ph/0508416v2.pdf))  
<a name="N8"></a>Nishikawa, H. & Kitamura, K. 2008, Journal of Computational Physics, 227, 2560 ([link](http://research.nianet.org/~hiro/My_papers/Nishikawa_Kitamura_JCP2008v227pp2560-2581Preprint.pdf))  

####Requirements
* Compiler support for C++11.
* [zlib](http://www.zlib.net): "A massively spiffy yet delicately unobtrusive compression library".  
* [Eigen](http://eigen.tuxfamily.org): "A template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms".  
* [Selene](https://github.com/jeremyong/Selene): "Simple C++11 friendly header-only bindings to Lua 5.2+".  
* [Lua5.2+](http://www.lua.org/): "A powerful, fast, lightweight, embeddable scripting language".  
