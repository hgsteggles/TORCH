/**
 * @file main.cpp
 */

#include "Torch.hpp"
#include "MPI_Wrapper.hpp"

#include <dirent.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include "DataPrinter.hpp"
#include "Logger.hpp"

#include "selene.h"

void parseParameters(const std::string& filename, TorchParameters& p);
void deleteFileContents(const std::string& myString);

int main (int argc, char** argv) {
	MPIW& mpihandler = MPIW::Instance(&argc, &argv);

	if (mpihandler.getRank() == 0)
		deleteFileContents("tmp/");
	try {
		TorchParameters tpars;
		mpihandler.serial([&] () { parseParameters( "refdata/parameters.lua", tpars ); });
		Torch torch;
		torch.initialise(tpars);
		torch.run();
	}
	catch (std::exception& e) {
		Logger<FileLogPolicy>::Instance().print<SeverityType::FATAL_ERROR>(e.what());
		MPIW::Instance().abort();
	}

	return 0;
}

std::string parseString(sel::Selector selector) {
	std::string str = selector;
	return str;
}

void parseParameters(const std::string& filename, TorchParameters& p) {
	Logger<FileLogPolicy>& logger = Logger<FileLogPolicy>::Instance();
	MPIW& mpihandler = MPIW::Instance();
	// Create new Lua state and load the lua libraries
	sel::State luaState{true};

	if (!luaState.Load(filename)) {
		throw std::runtime_error("ParseParameters: could not open lua file: " + filename);
	}
	else {
		p.nHI                  = (double)luaState["Parameters"]["Integration"]["ambient_no_density"];
		p.dscale               = (double)luaState["Parameters"]["Integration"]["density_scale"];
		p.pscale               = (double)luaState["Parameters"]["Integration"]["pressure_scale"];
		p.tscale               = (double)luaState["Parameters"]["Integration"]["time_scale"];
		p.spatialOrder         = (int)luaState["Parameters"]["Integration"]["spatial_order"];
		p.temporalOrder        = (int)luaState["Parameters"]["Integration"]["temporal_order"];
		p.tmax                 = (double)luaState["Parameters"]["Integration"]["simulation_time"];
		p.radiation_on         = (bool)luaState["Parameters"]["Integration"]["radiation_on"];
		p.cooling_on           = (bool)luaState["Parameters"]["Integration"]["cooling_on"];
		p.debug                = (bool)luaState["Parameters"]["Integration"]["debug"];
		p.outputDirectory      = parseString(luaState["Parameters"]["Integration"]["output_directory"]);
		p.initialConditions    = parseString(luaState["Parameters"]["Integration"]["initial_conditions"]);

		logger.print<SeverityType::DEBUG>("Data output file: ", p.outputDirectory);

		p.nd                   = (int)luaState["Parameters"]["Grid"]["no_dimensions"];
		p.ncells[0]            = (int)luaState["Parameters"]["Grid"]["no_cells_x"];
		p.ncells[1]            = (int)luaState["Parameters"]["Grid"]["no_cells_y"];
		p.ncells[2]            = (int)luaState["Parameters"]["Grid"]["no_cells_z"];
		p.sideLength           = (double)luaState["Parameters"]["Grid"]["side_length"];
		p.geometry             = parseString(luaState["Parameters"]["Grid"]["geometry"]);
		p.leftBC[0]            = parseString(luaState["Parameters"]["Grid"]["left_boundary_condition_x"]);
		p.leftBC[1]            = parseString(luaState["Parameters"]["Grid"]["left_boundary_condition_y"]);
		p.leftBC[2]            = parseString(luaState["Parameters"]["Grid"]["left_boundary_condition_z"]);
		p.rightBC[0]           = parseString(luaState["Parameters"]["Grid"]["right_boundary_condition_x"]);
		p.rightBC[1]           = parseString(luaState["Parameters"]["Grid"]["right_boundary_condition_y"]);
		p.rightBC[2]           = parseString(luaState["Parameters"]["Grid"]["right_boundary_condition_z"]);

		for (int i = 0; i < 3; ++i) {
			logger.print<SeverityType::DEBUG>("Boundary condition left[", i, "]: ", p.leftBC[i]);
			logger.print<SeverityType::DEBUG>("Boundary condition right[", i, "]: ", p.rightBC[i]);
		}

		p.heatCapacityRatio    = (double)luaState["Parameters"]["Hydrodynamics"]["gamma"];
		p.dfloor               = (double)luaState["Parameters"]["Hydrodynamics"]["density_floor"];
		p.pfloor               = (double)luaState["Parameters"]["Hydrodynamics"]["pressure_floor"];
		p.tfloor               = (double)luaState["Parameters"]["Hydrodynamics"]["temperature_floor"];
		p.riemannSolver        = parseString(luaState["Parameters"]["Hydrodynamics"]["riemann_solver"]);
		p.slopeLimiter         = parseString(luaState["Parameters"]["Hydrodynamics"]["slope_limiter"]);

		logger.print<SeverityType::DEBUG>("Using Riemann solver: ", p.riemannSolver);
		logger.print<SeverityType::DEBUG>("Using slope limiter: " + p.slopeLimiter);

		p.K1                   = (double)luaState["Parameters"]["Radiation"]["K1"];
		p.K2                   = (double)luaState["Parameters"]["Radiation"]["K2"];
		p.K3                   = (double)luaState["Parameters"]["Radiation"]["K3"];
		p.K4                   = (double)luaState["Parameters"]["Radiation"]["K4"];
		p.photoIonCrossSection = (double)luaState["Parameters"]["Radiation"]["photoion_cross_section"];
		p.alphaB               = (double)luaState["Parameters"]["Radiation"]["case_b_recombination_coeff"];
		p.tau0                 = (double)luaState["Parameters"]["Radiation"]["tau_0"];
		p.minX                 = (double)luaState["Parameters"]["Radiation"]["minimum_hii_fraction"];
		p.THI                  = (double)luaState["Parameters"]["Radiation"]["temperature_hi"];
		p.THII                 = (double)luaState["Parameters"]["Radiation"]["temperature_hii"];
		p.massFractionH        = (double)luaState["Parameters"]["Radiation"]["mass_fraction_hydrogen"];
		p.collisions_on        = (bool)luaState["Parameters"]["Radiation"]["collisions_on"];
		p.rt_coupling          = parseString(luaState["Parameters"]["Radiation"]["coupling"]);
		p.rt_scheme            = parseString(luaState["Parameters"]["Radiation"]["integration_scheme"]);

		logger.print<SeverityType::DEBUG>("Radiation/Hydrodynamics coupling: " + p.rt_coupling);
		logger.print<SeverityType::DEBUG>("Radiation integration scheme: " + p.rt_scheme);

		p.heatingAmplification = (double)luaState["Parameters"]["Thermodynamics"]["heating_amplification"];

		p.star_on              = (bool)luaState["Parameters"]["Star"]["on"];
		p.star_position[0]     = (int)luaState["Parameters"]["Star"]["cell_position_x"];
		p.star_position[1]     = (int)luaState["Parameters"]["Star"]["cell_position_y"];
		p.star_position[2]     = (int)luaState["Parameters"]["Star"]["cell_position_z"];
		p.faceSnap[0]          = (bool)luaState["Parameters"]["Star"]["snap_to_face_left_x"];
		p.faceSnap[1]          = (bool)luaState["Parameters"]["Star"]["snap_to_face_left_y"];
		p.faceSnap[2]          = (bool)luaState["Parameters"]["Star"]["snap_to_face_left_z"];
		p.photonEnergy         = (double)luaState["Parameters"]["Star"]["photon_energy"];
		p.photonRate           = (double)luaState["Parameters"]["Star"]["photon_rate"];
		p.windCellRadius       = (int)luaState["Parameters"]["Star"]["wind_radius_in_cells"];
		p.massLossRate         = (double)luaState["Parameters"]["Star"]["mass_loss_rate"];
		p.windVelocity         = (double)luaState["Parameters"]["Star"]["wind_velocity"];
		p.windTemperature      = (double)luaState["Parameters"]["Star"]["wind_temperature"];

		logger.print<SeverityType::DEBUG>("Star is on: " + p.star_on);
	}
}

////// deleteFileContents deletes all files within a directory
void deleteFileContents(const std::string& folder){
	struct dirent *next_file;
	DIR *dir; // These are data types defined in the "dirent" header.
	char filepath[256];
	dir = opendir(folder.c_str() );
	if (!dir) {
		throw std::runtime_error("deleteFileContents: directory " + folder + " does not exist.");
	}
	else {
		while ((next_file = readdir(dir) )) {
			// Build the full path for each file in the folder.
			sprintf(filepath, "%s/%s", folder.c_str(), next_file->d_name);
			remove(filepath);
		}
	}

	closedir(dir);
}
