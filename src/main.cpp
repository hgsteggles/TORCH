/**
 * @file main.cpp
 */

#include "Torch/Torch.hpp"
#include "MPI/MPI_Wrapper.hpp"
#include "IO/DataPrinter.hpp"
#include "IO/Logger.hpp"
#include "ParseLua.hpp"

#include "selene.h"

#include <dirent.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

void parseParameters(const std::string& filename, TorchParameters& p, int param_id);
void makeDirectory(const std::string& dir_name);
void deleteFileContents(const std::string& myString);
void copyConfigFile(const std::string& filename, const std::string& directory);

int main (int argc, char** argv) {
	MPIW& mpihandler = MPIW::Instance();

	try {
		TorchParameters tpars;
		std::string paramFile = "parameters.lua";
		int param_id = argc > 1 ? std::stoi(std::string(argv[1])) : 0;
		mpihandler.serial([&] () { parseParameters( paramFile, tpars, param_id ); });

		if (mpihandler.getRank() == 0) {
			makeDirectory(tpars.outputDirectory);
			deleteFileContents(tpars.outputDirectory);
			copyConfigFile("setup.lua", tpars.outputDirectory);
			copyConfigFile("parameters.lua", tpars.outputDirectory);
		}

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

void parseParameters(const std::string& filename, TorchParameters& p, int param_id) {
	std::string fullfilename = "refdata/" + filename;
	Logger<FileLogPolicy>& logger = Logger<FileLogPolicy>::Instance();
	// Create new Lua state and load the lua libraries
	sel::State luaState{true};

	if (!luaState.Load(fullfilename)) {
		throw std::runtime_error("ParseParameters: could not open lua file: " + fullfilename);
	}
	else {
		if (param_id != 0)
			luaState["param_set"](param_id);

		parseLuaVariable(luaState["Parameters"]["Integration"]["density_scale"], p.dscale);
		parseLuaVariable(luaState["Parameters"]["Integration"]["pressure_scale"], p.pscale);
		parseLuaVariable(luaState["Parameters"]["Integration"]["time_scale"], p.tscale);
		parseLuaVariable(luaState["Parameters"]["Integration"]["spatial_order"], p.spatialOrder);
		parseLuaVariable(luaState["Parameters"]["Integration"]["temporal_order"], p.temporalOrder);
		parseLuaVariable(luaState["Parameters"]["Integration"]["simulation_time"], p.tmax);
		parseLuaVariable(luaState["Parameters"]["Integration"]["radiation_on"], p.radiation_on);
		parseLuaVariable(luaState["Parameters"]["Integration"]["cooling_on"], p.cooling_on);
		parseLuaVariable(luaState["Parameters"]["Integration"]["debug"], p.debug);
		parseLuaVariable(luaState["Parameters"]["Integration"]["setup_id"], p.setupID);
		parseLuaVariable(luaState["Parameters"]["Integration"]["output_directory"], p.outputDirectory);
		parseLuaVariable(luaState["Parameters"]["Integration"]["initial_conditions"], p.initialConditions);

		parseLuaVariable(luaState["Parameters"]["Grid"]["no_dimensions"], p.nd);
		parseLuaVariable(luaState["Parameters"]["Grid"]["no_cells_x"], p.ncells[0]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["no_cells_y"], p.ncells[1]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["no_cells_z"], p.ncells[2]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["side_length"], p.sideLength);
		parseLuaVariable(luaState["Parameters"]["Grid"]["geometry"], p.geometry);
		parseLuaVariable(luaState["Parameters"]["Grid"]["left_boundary_condition_x"], p.leftBC[0]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["left_boundary_condition_y"], p.leftBC[1]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["left_boundary_condition_z"], p.leftBC[2]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["right_boundary_condition_x"], p.rightBC[0]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["right_boundary_condition_y"], p.rightBC[1]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["right_boundary_condition_z"], p.rightBC[2]);

		parseLuaVariable(luaState["Parameters"]["Grid"]["Patch"]["filename"], p.patchfilename);
		parseLuaVariable(luaState["Parameters"]["Grid"]["Patch"]["offset_x"], p.patchoffset[0]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["Patch"]["offset_y"], p.patchoffset[1]);
		parseLuaVariable(luaState["Parameters"]["Grid"]["Patch"]["offset_z"], p.patchoffset[2]);

		parseLuaVariable(luaState["Parameters"]["Hydrodynamics"]["gamma"], p.heatCapacityRatio);
		parseLuaVariable(luaState["Parameters"]["Hydrodynamics"]["density_floor"], p.dfloor);
		parseLuaVariable(luaState["Parameters"]["Hydrodynamics"]["pressure_floor"], p.pfloor);
		parseLuaVariable(luaState["Parameters"]["Hydrodynamics"]["temperature_floor"], p.tfloor);
		parseLuaVariable(luaState["Parameters"]["Hydrodynamics"]["riemann_solver"], p.riemannSolver);
		parseLuaVariable(luaState["Parameters"]["Hydrodynamics"]["slope_limiter"], p.slopeLimiter);

		parseLuaVariable(luaState["Parameters"]["Radiation"]["K1"], p.K1);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["K2"], p.K2);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["K3"], p.K3);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["K4"], p.K4);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["photoion_cross_section"], p.photoIonCrossSection);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["case_b_recombination_coeff"], p.alphaB);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["tau_0"], p.tau0);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["minimum_hii_fraction"], p.minX);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["temperature_hi"], p.THI);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["temperature_hii"], p.THII);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["mass_fraction_hydrogen"], p.massFractionH);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["collisions_on"], p.collisions_on);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["coupling"], p.rt_coupling);
		parseLuaVariable(luaState["Parameters"]["Radiation"]["integration_scheme"], p.rt_scheme);

		parseLuaVariable(luaState["Parameters"]["Thermodynamics"]["thermo_hii_switch"], p.thermoHII_Switch);
		parseLuaVariable(luaState["Parameters"]["Thermodynamics"]["heating_amplification"], p.heatingAmplification);
		parseLuaVariable(luaState["Parameters"]["Thermodynamics"]["thermo_subcycling"], p.thermoSubcycling);
		parseLuaVariable(luaState["Parameters"]["Thermodynamics"]["min_temp_initial_state"], p.minTempInitialState);

		parseLuaVariable(luaState["Parameters"]["Star"]["on"], p.star_on);
		parseLuaVariable(luaState["Parameters"]["Star"]["cell_position_x"], p.star_position[0]);
		parseLuaVariable(luaState["Parameters"]["Star"]["cell_position_y"], p.star_position[1]);
		parseLuaVariable(luaState["Parameters"]["Star"]["cell_position_z"], p.star_position[2]);
		parseLuaVariable(luaState["Parameters"]["Star"]["snap_to_face_left_x"], p.faceSnap[0]);
		parseLuaVariable(luaState["Parameters"]["Star"]["snap_to_face_left_y"], p.faceSnap[1]);
		parseLuaVariable(luaState["Parameters"]["Star"]["snap_to_face_left_z"], p.faceSnap[2]);
		parseLuaVariable(luaState["Parameters"]["Star"]["photon_energy"], p.photonEnergy);
		parseLuaVariable(luaState["Parameters"]["Star"]["photon_rate"], p.photonRate);
		parseLuaVariable(luaState["Parameters"]["Star"]["wind_radius_in_cells"], p.windCellRadius);
		parseLuaVariable(luaState["Parameters"]["Star"]["mass_loss_rate"], p.massLossRate);
		parseLuaVariable(luaState["Parameters"]["Star"]["wind_velocity"], p.windVelocity);
		parseLuaVariable(luaState["Parameters"]["Star"]["wind_temperature"], p.windTemperature);

		logger.print<SeverityType::DEBUG>("Star is on: ", p.star_on);
	}
}

void makeDirectory(const std::string& dir_name) {
	struct stat st = {0};
	if (stat(dir_name.c_str(), &st) == -1)
		mkdir(dir_name.c_str(), 0755);
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

void copyConfigFile(const std::string& filename, const std::string& directory) {
	std::string configDirectory = "refdata/";
	std::ifstream  src(configDirectory + filename, std::ios::binary);
	std::ofstream  dst(directory + filename, std::ios::binary);
	dst << src.rdbuf();
}
