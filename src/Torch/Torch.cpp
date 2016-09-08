#include "Torch.hpp"
#include "IO/ProgressBar.hpp"
#include "Fluid/GridCell.hpp"
#include "Constants.hpp"
#include "IO/Logger.hpp"
#include "Misc/Timer.hpp"
#include "IO/DataReader.hpp"

#include <chrono>
#include <fstream>
#include <string>
#include <iostream>
#include <assert.h>
#include <cmath>

#include "selene/include/selene.h"

int stepIDFromFilename(const std::string& filename) {
	std::size_t lastindex = filename.find_last_of(".");
	std::string rawname = (lastindex != std::string::npos) ? filename.substr(0, lastindex) : filename;

	std::size_t underscore_index = rawname.find_last_of("_");
	std::string stepno_string = (underscore_index != std::string::npos) ? rawname.substr(underscore_index+1) : "-1";

	return std::stoi(stepno_string);
}

void Torch::initialise(TorchParameters p) {
	consts = std::make_shared<Constants>();

	// Initialise the scalings (scaling physical units to code units (to reduce chance of arithmetic underflow/overflow). 
	consts->initialise_DPT(p.dscale, p.pscale, p.tscale);

	// Initialise code parameters.
	p.initialise(consts);

	// Read grid geometry from initial conditions data file if one exists to set up initial grid data structure.
	DataParameters datap;
	if (p.initialConditions.compare("") != 0) {
		datap = DataReader::readDataParameters(p.initialConditions);
		p.ncells = datap.ncells;
		p.sideLength = consts->converter.toCodeUnits(datap.sideLength, 0, 1, 0);
		p.nd = datap.nd;
	}

	// Forward parameters to Constants object.
	consts->nd = p.nd;
	consts->dfloor = p.dfloor;
	consts->pfloor = p.pfloor;
	consts->tfloor = p.tfloor;

	// Initialise IO with output directory and consts (which includes unit conversion info).
	inputOutput.initialise(consts, p.outputDirectory);

	// Set up grid data structure using geometry info read in earlier.
	fluid.initialise(consts, p.getFluidParameters());
	fluid.initialiseGrid(p.getGridParameters(), p.getStarParameters());
	fluid.getGrid().currentTime = consts->converter.toCodeUnits(datap.time, 0, 0, 1);

	// Forward hydrodynamics parameters.
	hydrodynamics.initialise(consts);

	// Try to set up RiemannSolver and SlopeLimiter with strings passed in parameters.lua - if invalid the default is used and a warning is issued to the log file.
	try {
		hydrodynamics.setRiemannSolver(std::move(RiemannSolverFactory::create(p.riemannSolver, p.nd)));
	}
	catch (std::exception& e) {
		Logger<FileLogPolicy>::Instance().print<SeverityType::WARNING>(e.what());
	}

	try {
		hydrodynamics.setSlopeLimiter(std::move(SlopeLimiterFactory::create(p.slopeLimiter)));
	}
	catch (std::exception& e) {
		Logger<FileLogPolicy>::Instance().print<SeverityType::WARNING>(e.what());
	}

	// Forward parameters to Radiation object.
	radiation.initialise(consts, p.getRadiationParameters());
	// Forward parameters to Thermodynamics object.
	thermodynamics.initialise(consts, p.getThermoParameters());

	// Forward parameters to this object.
	initialConditions = p.initialConditions;
	radiation_on = p.radiation_on;
	cooling_on = p.cooling_on;
	debug = p.debug;
	spatialOrder = p.spatialOrder;
	temporalOrder = p.temporalOrder;
	tmax = p.tmax;
	dt_max = p.dt_max;
	dfloor = p.dfloor;
	pfloor = p.pfloor;
	tfloor = p.tfloor;

	steps = 0;
	stepCounter = 0;

	if (initialConditions.compare("") != 0) {
		std::cout << initialConditions << std::endl;
		//setUp(initialConditions);
		DataReader::readGrid(p.initialConditions, datap, fluid);
		Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>("Torch::initialise: Grid read from file.");
		stepstart = stepIDFromFilename(p.initialConditions);
	}
	else {
		// Set up initial grid state using the setup.lua file.
		setUpLua(p.setupFile);
	}
	if (p.patchfilename.compare("") != 0)
		DataReader::patchGrid(p.patchfilename, p.patchoffset, fluid);

	// Initialise the minimum temperature of cells given the initial temperature field if this is turned on in the parameters.lua file.
	thermodynamics.initialiseMinTempField(fluid);

	// Convert cell data to code units, fix any broken primitive variables and calculate the conservative variables.
	toCodeUnits();
	fluid.fixPrimitives();
	fluid.globalUfromQ();

	// Initialise the path lengths, shell volumes, and nearest neighbour weights for use with the radiative transfer.
	radiation.initField(fluid);

	// Warn the user if the reverse shock of the star is within or close to the injection radius.
	if (p.star_on && p.windCellRadius > 0) {
		Star& star = fluid.getStar();
		if (star.core == Star::Location::HERE) {
			Grid& grid = fluid.getGrid();

			double edot = 0.5 * star.massLossRate * star.windVelocity * star.windVelocity;
			double pre = grid.getCell(grid.locate((int)star.xc[0], (int)star.xc[1], (int)star.xc[2])).Q[UID::PRE];
			double reverse2 = std::sqrt(2.0 * edot * fluid.getStar().massLossRate) / (4.0 * consts->pi * pre);
			double reverse = std::sqrt(reverse2) / fluid.getGrid().dx[0];
			if (reverse < 5 + p.windCellRadius) {
				std::cout << "Warning: reverse shock within or close to wind injection region:" << std::endl;
				std::cout << "         [rs = " << reverse << ", wir = " << p.windCellRadius << "]" << std::endl;
			}
		}
	}

	Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>("Torch::initialise: initial setup complete.");
}

void Torch::toCodeUnits() {
	for (GridCell& cell : fluid.getGrid().getIterable("GridCells")) {
		cell.Q[UID::DEN] = consts->converter.toCodeUnits(cell.Q[UID::DEN], 1, -3, 0);
		cell.Q[UID::PRE] = consts->converter.toCodeUnits(cell.Q[UID::PRE], 1, -1, -2);
		for (int idim = 0; idim < consts->nd; ++idim)
			cell.Q[UID::VEL+idim] = consts->converter.toCodeUnits(cell.Q[UID::VEL+idim], 0, 1, -1);
		for (int idim = 0; idim < consts->nd; ++idim)
			cell.GRAV[idim] = consts->converter.toCodeUnits(cell.GRAV[idim], 1, -2, -2);
	}
}

void Torch::setUp(std::string filename) {
	MPIW& mpihandler = MPIW::Instance();
	double ignore;
	mpihandler.serial([&] () {
		std::fstream myfile(filename, std::ios_base::in);
		myfile >> fluid.getGrid().currentTime >> ignore >> ignore >> ignore;
		fluid.getGrid().currentTime = consts->converter.toCodeUnits(fluid.getGrid().currentTime, 0, 0, 1);

		int skip = mpihandler.getRank()*fluid.getGrid().ncells[0]*fluid.getGrid().ncells[1]*fluid.getGrid().ncells[2]/mpihandler.nProcessors();

		for (int i = 0; i < skip; ++i) {
			for (int idim = 0; idim < consts->nd; ++idim)
				myfile >> ignore >> ignore;
			myfile >> ignore >> ignore >> ignore;
		}

		for (GridCell& cell : fluid.getGrid().getIterable("GridCells")) {
			for (int idim = 0; idim < consts->nd; ++idim)
				myfile >> ignore;

			myfile >> cell.Q[UID::DEN];
			myfile >> cell.Q[UID::PRE];
			myfile >> cell.Q[UID::HII];

			for (int idim = 0; idim < consts->nd; ++idim) {
				myfile >> cell.Q[UID::VEL+idim];
			}
			cell.heatCapacityRatio = fluid.heatCapacityRatio;
		}
		myfile.close();
	});

	Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>("Torch::setUp(initialConditionsFile) complete.");
}

void Torch::setUpLua(std::string filename) {
	MPIW& mpihandler = MPIW::Instance();
	Grid& grid = fluid.getGrid();

	mpihandler.print("Reading lua config file: " + filename);

	mpihandler.serial([&] () {
		// Create new Lua state and load the lua libraries
		sel::State luaState{true};
		bool hasLoaded = luaState.Load(filename);

		if (!hasLoaded)
			std::cout << "SetUpLua: could not open lua file: " << filename << std::endl;
		else {
			for (GridCell& cell : grid.getIterable("GridCells")) {
				std::array<double, 3> xc, xs;
				for (int i = 0; i < 3; ++i) {
					xc[i] = consts->converter.fromCodeUnits(cell.xc[i]*grid.dx[i], 0, 1, 0);
					xs[i] = consts->converter.fromCodeUnits(fluid.getStar().xc[i]*grid.dx[i], 0, 1, 0);
				}

				sel::tie(cell.Q[UID::DEN],
						cell.Q[UID::PRE],
						cell.Q[UID::HII],
						cell.Q[UID::VEL],
						cell.Q[UID::VEL+1],
						cell.Q[UID::VEL+2],
						cell.GRAV[0],
						cell.GRAV[1],
						cell.GRAV[2])
					= luaState["initialise"](xc[0], xc[1], xc[2], xs[0], xs[1], xs[2]);

				cell.heatCapacityRatio = fluid.heatCapacityRatio;
			}
		}
	});
}

void Torch::run() {
	MPIW& mpihandler = MPIW::Instance();
	Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>("Torch::run() initial conditions set up.");

	double initTime = fluid.getGrid().currentTime;
	ProgressBar progBar = ProgressBar(tmax, 1, "Marching solution", debug);

	fluid.globalQfromU();
	fluid.fixPrimitives();

	progBar.update(initTime, dt_max, mpihandler.getRank() == 0);
	inputOutput.print2D(std::to_string((int)(100.0*initTime/tmax + 0.5)), initTime, fluid.getGrid());

	Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>("Torch::run() first data dump complete.");

	activeComponents.push_back(ComponentID::HYDRO);
	if (cooling_on)
		activeComponents.push_back(ComponentID::THERMO);
	if (radiation_on)
		activeComponents.push_back(ComponentID::RAD);

	Timer timer;
	timer.start();

	bool isFinalPrintOn = true;

	while (fluid.getGrid().currentTime < tmax && !m_isQuitting) {
		// Find the time until the next data snapshot. Print if it has passed.
		double dt_nextCheckpoint = dt_max;
		bool print_now = progBar.update(fluid.getGrid().currentTime, dt_nextCheckpoint, mpihandler.getRank() == 0);
		if (print_now) {
			int step = (int)(100.0*fluid.getGrid().currentTime/tmax + 0.5);
			inputOutput.print2D(std::to_string(step), fluid.getGrid().currentTime, fluid.getGrid());
			isFinalPrintOn = (step != 100);
		}

		// Perform full integration time-step of all physics sub-problems.
		fluid.getGrid().deltatime = fullStep(dt_nextCheckpoint);
		fluid.getGrid().currentTime += fluid.getGrid().deltatime;
		++steps;
	}
	progBar.end(mpihandler.getRank() == 0);

	if (isFinalPrintOn) {
		inputOutput.print2D(std::to_string(100), fluid.getGrid().currentTime, fluid.getGrid());
	}

	mpihandler.barrier();
	if (mpihandler.getRank() == 0)
		std::cout << "MARCH: Took " << timer.formatTime(timer.getTicks()) << '\n';
}

double Torch::calculateTimeStep() {
	static bool first_time = true;
	double dt;
	if (first_time) {
		dt = dt_max*1.0e-20;
		first_time = false;
	}
	else {
		double dt_hydro = hydrodynamics.calculateTimeStep(dt_max, fluid);
		double dt_rad = dt_hydro;
		double dt_thermo = dt_hydro;
		if (radiation_on)
			dt_rad = radiation.calculateTimeStep(dt_max, fluid);
		if (cooling_on)
			dt_thermo = thermodynamics.calculateTimeStep(dt_max, fluid);
		dt = std::min(std::min(dt_hydro, dt_rad), dt_thermo);

		if (debug) {
			double thyd = 100.0*dt_hydro/tmax;
			thyd = MPIW::Instance().minimum(thyd);
			double trad = 100.0*dt_rad/tmax;
			trad = MPIW::Instance().minimum(trad);
			double ttherm = 100.0*dt_thermo/tmax;
			ttherm = MPIW::Instance().minimum(ttherm);
			if (MPIW::Instance().getRank() == 0)
				std::cout << "thyd = " << thyd << ", trad = " << trad << ", ttherm = " << ttherm << '\n';
			if (thyd <= 1.0e-6 || thyd <= 1.0e-6 || thyd <= 1.0e-6)
				m_isQuitting = true;
		}
	}
	dt = MPIW::Instance().minimum(dt);
	inputOutput.reduceToPrint(fluid.getGrid().currentTime, dt);
	fluid.getGrid().deltatime = dt;
	return dt;
}

Integrator& Torch::getComponent(ComponentID id) {
	switch(id) {
	case ComponentID::RAD:
		return radiation;
	case ComponentID::THERMO:
		return thermodynamics;
	default:
		return hydrodynamics;
	}
}

void Torch::subStep(double dt, bool hasCalculatedHeatFlux, Integrator& comp) {
	checkValues(comp.getComponentName() + "before");
	if (!hasCalculatedHeatFlux) {
		fluid.globalQfromU();
		fluid.fixPrimitives();
		comp.preTimeStepCalculations(fluid);
	}
	comp.integrate(dt, fluid);
	comp.updateSourceTerms(dt, fluid);
	fluid.advSolution(dt);
	fluid.fixSolution();
	checkValues(comp.getComponentName() + " after");
}

void Torch::hydroStep(double dt, bool hasCalculatedHeatFlux) {
	checkValues("hydro before");
	fluid.globalWfromU();
	if (!hasCalculatedHeatFlux) {
		fluid.globalQfromU();
		fluid.fixPrimitives();
		hydrodynamics.preTimeStepCalculations(fluid);
	}
	hydrodynamics.integrate(dt, fluid);
	hydrodynamics.updateSourceTerms(dt, fluid);

	fluid.advSolution(dt/2.0);
	fluid.fixSolution();

	// Corrector.
	fluid.globalQfromU();
	fluid.globalUfromW();
	hydrodynamics.integrate(dt, fluid);
	hydrodynamics.updateSourceTerms(dt, fluid);
	fluid.advSolution(dt);
	fluid.fixSolution();
}

double Torch::fullStep(double dt_nextCheckPoint) {
	fluid.globalQfromU();
	fluid.fixPrimitives();
	if (cooling_on)
		thermodynamics.preTimeStepCalculations(fluid);
	if (radiation_on)
		radiation.preTimeStepCalculations(fluid);

	double dt = std::min(dt_nextCheckPoint, calculateTimeStep());

	int ncomps = activeComponents.size();

	if (ncomps == 1) {
		hydroStep(dt, true);
		return dt;
	}

	stepCounter = (stepCounter+1)%ncomps;

	for (int i = 0; i < ncomps; ++i) {
		double h = (i == ncomps-1) ? 1.0 : 0.5;
		subStep(h*dt, i == 0, getComponent(activeComponents[(i+stepCounter)%ncomps]));
	}

	for (int i = ncomps-2; i >= 0; --i) {
		subStep(dt/2.0, false, getComponent(activeComponents[(i+stepCounter)%ncomps]));
	}

	return dt;
}

void Torch::checkValues(std::string componentname) {
	bool error = false;
	for (GridCell& cell : fluid.getGrid().getIterable("GridCells")) {
		for (int i = 0; i < UID::N; ++i) {
			if (cell.U[i] != cell.U[i] || std::isinf(cell.U[i]) || cell.Q[UID::DEN] == 0 || cell.Q[UID::PRE] == 0) {
				error = true;
				break;
			}
		}
		if (error)
			break;
	}
	if (error) {
		for (GridCell& cell : fluid.getGrid().getIterable("GridCells")) {
			if (std::abs(cell.Q[UID::VEL+0]) > 1e50 || std::abs(cell.Q[UID::VEL+1]) > 1e50) {
				std::cout << '\n' << componentname << " produced an error.\n";
				cell.printInfo();
			}
		}
		exit(0);
	}
}

/*
void Torch::checkValues(std::string componentname) {
	for (GridCell& cell : fluid.getGrid().getCausalCells()) {
		bool error = false;
		for (int i = 0; i < UID::N; ++i) {
			if (cell.U[i] != cell.U[i] || std::isinf(cell.U[i])) {
				error = true;
				break;
			}
		}
		if (error) {
			std::cout << '\n' << componentname << " produced an error.\n";
			if (componentname.compare("Hydrodynamics") == 0) {
				std::cout << "Possible Cause: is the RiemannSolver stable?\n";
				std::cout << "Possible Cause: is the SlopeLimiter TVD?\n";
				std::cout << "Possible Cause: are the boundaries correctly linked to grid? \n";
			}
			else if (componentname.compare("Radiation") == 0 || componentname.compare("Thermodynamics") == 0) {
				std::cout << "Possible Cause: are the column densities calculated properly?\n";
			}
			std::cout << "Possible Cause: have the variables been appropriately floored?\n";
			cell.printInfo();
			exit(EXIT_FAILURE);
		}
	}
}
 */
