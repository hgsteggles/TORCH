/**
 * @file main.cpp
 */

#include "Torch/Torch.hpp"
#include "MPI/MPI_Wrapper.hpp"
#include "IO/DataPrinter.hpp"
#include "IO/Logger.hpp"
#include "IO/ParseLua.hpp"
#include "IO/FileManagement.hpp"

#include "selene/include/selene.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>

std::string readOutputDirectory(const std::string& paramfilename);
void parseParameters(const std::string& filename, TorchParameters& p);
void showUsage();

int main (int argc, char** argv) {
	MPIW& mpihandler = MPIW::Instance();

	std::string paramFile = "config/torch-config.lua";
	std::string setupFile = "config/torch-setup.lua";

	// Parse parameters
	if (argc > 3) {
		showUsage();
		exit(0);
	}
	if (argc > 1) {
		for (int iarg = 1; iarg < argc; ++iarg) {
			std::string arg = argv[iarg];
			std::string prefix1("--paramfile=");
			std::string prefix2("--setupfile=");

			if (!arg.compare(0, prefix1.size(), prefix1))
				paramFile = arg.substr(prefix1.size()).c_str();
			else if (!arg.compare(0, prefix2.size(), prefix2))
				setupFile = arg.substr(prefix2.size()).c_str();
			else {
				showUsage();
				exit(0);
			}
		}
	}

	TorchParameters tpars;

	try {
		mpihandler.serial([&] () {
			tpars.outputDirectory = readOutputDirectory(paramFile);
		});

		if (mpihandler.getRank() == 0) {
			FileManagement::makeDirectoryPath(tpars.outputDirectory);
			FileManagement::deleteFileContents(tpars.outputDirectory);
			FileManagement::makeDirectoryPath(tpars.outputDirectory + "/log");
			FileManagement::deleteFileContents(tpars.outputDirectory + "/log");
			FileManagement::copyConfigFile(setupFile, tpars.outputDirectory);
			FileManagement::copyConfigFile(paramFile, tpars.outputDirectory);
		}
	}
	catch (std::exception& e) {
		std::cout << e.what() << std::endl;
		MPIW::Instance().abort();
	}

	Logger<FileLogPolicy>& logger = Logger<FileLogPolicy>::Instance(tpars.outputDirectory);

	try {
		mpihandler.serial([&] () {
			parseParameters( paramFile, tpars);
			tpars.setupFile = setupFile;
		});

		Torch torch;
		torch.initialise(tpars);
		torch.run();
	}
	catch (std::exception& e) {
		logger.print<SeverityType::FATAL_ERROR>(e.what());
		std::cout << "Program encountered a fatal error. See ./" << tpars.outputDirectory << "/log/torch.log*" << " for more details." << std::endl;
		MPIW::Instance().abort();
	}

	return 0;
}

std::string readOutputDirectory(const std::string& paramfilename) {
	std::string fullfilename = paramfilename;
	std::string outputDir = "";
	// Create new Lua state and load the lua libraries
	sel::State luaState{true};

	if (!luaState.Load(fullfilename)) {
		throw std::runtime_error("ParseParameters: could not open lua file: " + fullfilename);
	}
	else {
		parseLuaVariable(luaState["Parameters"]["Integration"]["output_directory"], outputDir);
	}

	return outputDir;
}

void parseParameters(const std::string& filename, TorchParameters& p) {
	std::string fullfilename = filename;
	Logger<FileLogPolicy>& logger = Logger<FileLogPolicy>::Instance();
	// Create new Lua state and load the lua libraries
	sel::State luaState{true};

	if (!luaState.Load(fullfilename)) {
		throw std::runtime_error("ParseParameters: could not open lua file: " + fullfilename);
	}
	else {
		parseLuaVariable(luaState["Parameters"]["Integration"]["density_scale"], p.dscale);
		parseLuaVariable(luaState["Parameters"]["Integration"]["pressure_scale"], p.pscale);
		parseLuaVariable(luaState["Parameters"]["Integration"]["time_scale"], p.tscale);
		parseLuaVariable(luaState["Parameters"]["Integration"]["spatial_order"], p.spatialOrder);
		parseLuaVariable(luaState["Parameters"]["Integration"]["temporal_order"], p.temporalOrder);
		parseLuaVariable(luaState["Parameters"]["Integration"]["simulation_time"], p.tmax);
		parseLuaVariable(luaState["Parameters"]["Integration"]["radiation_on"], p.radiation_on);
		parseLuaVariable(luaState["Parameters"]["Integration"]["cooling_on"], p.cooling_on);
		parseLuaVariable(luaState["Parameters"]["Integration"]["debug"], p.debug);
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

void showUsage() {
	std::cout << "torch [--paramfile=<filename>] [--setupfile=<filename>]" << std::endl;
}
