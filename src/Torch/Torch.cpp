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

#include "selene.h"

void Torch::initialise(TorchParameters p) {
	MPIW& mpihandler = MPIW::Instance();

	consts = std::make_shared<Constants>();
	consts->initialise_DPT(p.dscale, p.pscale, p.tscale);

	p.initialise(consts);

	/*
	mpihandler.serial([&] () {
		if (p.initialConditions.compare("") != 0) {
			std::ifstream myfile(p.initialConditions, std::ios_base::in);
			if (!myfile)
				throw std::runtime_error("Torch::initialise: invalid input file " + p.initialConditions + ".");
			double t;
			int a, b, c;
			myfile >> t >> a >> b >> c;
			p.ncells[0] = a;
			p.ncells[1] = b;
			p.ncells[2] = c;
			myfile.close();
		}
	});

	if (p.ncells[1] == 1)
		p.nd = 1;
	else if (p.ncells[2] == 1)
		p.nd = 2;
	else
		p.nd = 3;
	*/

	DataParameters datap;
	if (p.initialConditions.compare("") != 0) {
		datap = DataReader::readDataParameters(p.initialConditions);
		p.ncells = datap.ncells;
		p.sideLength = datap.sideLength;
		p.nd = datap.nd;
	}

	consts->nd = p.nd;

	consts->dfloor = p.dfloor;
	consts->pfloor = p.pfloor;
	consts->tfloor = p.tfloor;

	inputOutput.initialise(consts, p.outputDirectory);

	fluid.initialise(consts, p.getFluidParameters());
	fluid.initialiseGrid(p.getGridParameters(), p.getStarParameters());
	fluid.getGrid().currentTime = consts->converter.toCodeUnits(datap.time, 0, 0, 1);

	hydrodynamics.initialise(consts);

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

	radiation.initialise(consts, p.getRadiationParameters());

	thermodynamics.initialise(consts, p.getThermoParameters());

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

	if (initialConditions.compare("") == 0)
		setUpLua("refdata/setup.lua", p.setupID);
	else {
		//setUp(initialConditions);
		DataReader::readGrid(p.initialConditions, datap, fluid);
		Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>("Torch::setUp(initialConditionsFile) complete.");
	}
	if (p.patchfilename.compare("") != 0)
		DataReader::patchGrid(p.patchfilename, p.patchoffset, fluid);
	toCodeUnits();
	fluid.fixPrimitives();
	fluid.globalUfromQ();

	radiation.initField(fluid);
}

void Torch::toCodeUnits() {
	for (GridCell& cell : fluid.getGrid().getCells()) {
		cell.Q[UID::DEN] = consts->converter.toCodeUnits(cell.Q[UID::DEN], 1, -3, 0);
		cell.Q[UID::PRE] = consts->converter.toCodeUnits(cell.Q[UID::PRE], 1, -1, -2);
		for (int idim = 0; idim < consts->nd; ++idim)
			cell.Q[UID::VEL+idim] = consts->converter.toCodeUnits(cell.Q[UID::VEL+idim], 0, 1, -1);
		for (int idim = 0; idim < consts->nd; ++idim)
			cell.GRAV[idim] = consts->converter.toCodeUnits(cell.GRAV[idim], 0, 1, -2);
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

		for (GridCell& cell : fluid.getGrid().getCells()) {
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

void Torch::setUpLua(std::string filename, int setupID) {
	MPIW& mpihandler = MPIW::Instance();

	mpihandler.print("Reading lua config file: " + filename);

	mpihandler.serial([&] () {
		// Create new Lua state and load the lua libraries
		sel::State luaState{true};
		bool hasLoaded = luaState.Load(filename);

		if (!hasLoaded)
			std::cout << "SetUpLua: could not open lua file: " << filename << std::endl;
		else {
			if (setupID != 0)
				luaState["setup_set"](setupID);

			for (GridCell& cell : fluid.getGrid().getCells()) {
				std::array<double, 3> xc, xs;
				for (int i = 0; i < 3; ++i) {
					xc[i] = consts->converter.fromCodeUnits(cell.xc[i]*fluid.getGrid().dx[i], 0, 1, 0);
					xs[i] = consts->converter.fromCodeUnits(fluid.getStar().xc[i]*fluid.getGrid().dx[i], 0, 1, 0);
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
	Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>("Torch::setUpLua() complete.");
}

void Torch::run() {
	MPIW& mpihandler = MPIW::Instance();
	Logger<FileLogPolicy>::Instance().print<SeverityType::DEBUG>("Torch::run() initial conditions set up.");

	double initTime = fluid.getGrid().currentTime;
	ProgressBar progBar = ProgressBar(tmax-initTime, 1, "Marching solution", debug);

	fluid.globalQfromU();
	fluid.fixPrimitives();

	//hydrodynamics.fixIC(fluid);

	inputOutput.print2D(std::to_string(0), fluid.getGrid().currentTime, fluid.getGrid());
	inputOutput.printSTARBENCH(radiation, hydrodynamics, fluid);

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
		double dt_nextCheckpoint = dt_max;
		bool print_now = progBar.update(fluid.getGrid().currentTime-initTime, dt_nextCheckpoint, mpihandler.getRank() == 0);
		if (print_now) {
			int step = (int)(100.0*(fluid.getGrid().currentTime-initTime)/(tmax-initTime) + 0.5);
			inputOutput.print2D(std::to_string(step), fluid.getGrid().currentTime, fluid.getGrid());
			thermodynamics.fillHeatingArrays(fluid);
			inputOutput.printVariables(step, fluid.getGrid().currentTime, fluid.getGrid());
			isFinalPrintOn = (step != 100);
		}

		fluid.getGrid().deltatime = fullStep(dt_nextCheckpoint);
		fluid.getGrid().currentTime += fluid.getGrid().deltatime;
		++steps;
	}
	progBar.end(mpihandler.getRank() == 0);

	if (isFinalPrintOn) {
		thermodynamics.fillHeatingArrays(fluid);
		inputOutput.printVariables(100, fluid.getGrid().currentTime, fluid.getGrid());
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
	for (GridCell& cell : fluid.getGrid().getCells()) {
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
		for (GridCell& cell : fluid.getGrid().getCells()) {
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
