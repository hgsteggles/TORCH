/** Provides the Torch class.
 *
 * @file Torch.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef TORCH_HPP_
#define TORCH_HPP_

#include <memory>
#include <string>
#include <vector>

#include "Fluid/Fluid.hpp"
#include "Integrators/Hydro.hpp"
#include "Integrators/Radiation.hpp"
#include "Integrators/Riemann.hpp"
#include "Integrators/SlopeLimiter.hpp"
#include "Integrators/Thermodynamics.hpp"
#include "IO/DataPrinter.hpp"
#include "Parameters.hpp"

//#include "Star.hpp"

class Constants;

enum class ComponentID : unsigned int {HYDRO, RAD, THERMO};

/**
 * @class Torch
 *
 * @brief Class that sets up a problem and integrates the solution to the state of the Fluid.
 *
 * @version 0.8, 24/11/2014
 */
class Torch {
public:
	void initialise(TorchParameters tparams);
	void run();

	void setRiemannSolver(std::unique_ptr<RiemannSolver> riemannSolver);
	void setSlopeLimiter(std::unique_ptr<SlopeLimiter> slopeLimiter);

	bool isQuitting();
	TorchParameters getRemapParameters();
private:
	std::shared_ptr<Constants> consts = nullptr;
	DataPrinter inputOutput; //!< Module for input/output.
	Fluid fluid;
	Hydrodynamics hydrodynamics;
	Radiation radiation;
	Thermodynamics thermodynamics;

	TorchParameters remapParameters;
	std::string initialConditions = "";
	std::vector<ComponentID> activeComponents;
	bool radiation_on = false;
	bool cooling_on = false;
	bool debug = false;
	unsigned int spatialOrder = 0;
	unsigned int temporalOrder = 0;
	double tmax = 0;
	double dt_max = 0;
	double dfloor = 0;
	double pfloor = 0;
	double tfloor = 0;
	long steps = 0;
	long stepstart = 0;
	int stepCounter = 0;
	int m_customPrintID = 0;

	bool m_isQuitting = false;

	void toCodeUnits();
	void setUp(std::string filename);
	void setUpLua(std::string filename, int setupID);
	double calculateTimeStep();
	Integrator& getComponent(ComponentID id);
	void hydroStep(double dt, bool hasCalculatedHeatFlux);
	void subStep(double dt, bool hasCalculatedHeatFlux, Integrator& comp);
	double fullStep(double dt_nextCheckPoint);
	void checkValues(std::string componentname);
};

#endif // TORCH_HPP_
