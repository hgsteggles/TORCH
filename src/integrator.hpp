/**
 * @file integrator.hpp
 *
 * @author Harrison Steggles
 * @date 04/02/2014 - initial verision. Encapsulating program into integrator class.
 * @date 12/02/2014 - command line progress bar added to march method. Removed previous way to display progress.
 * @date 01/05/2014 - Integrator now has an IntegrationParameters member.
 * @date 28/05/2014 - march steps the fluid hydro/radiation variables using Strang Splitting instead of weird
 * method I was implementing before.
 */

#ifndef INTEGRATOR_HPP_
#define INTEGRATOR_HPP_

#include <string>

class Grid3D;
class Radiation;
class HydroDynamics;
class InputOutput;
class MPIHandler;
class IntegrationParameters;
class GridParameters;
class RadiationParameters;
class HydroParameters;
class PrintParameters;
class Scalings;

class Integrator {
public:
	IntegrationParameters& iparams; //!< Contains integration parameters.
	Grid3D* grid; //!< The grid structure that holds integration cells.
	Radiation* rad; //!< Module for radiative transfer.
	HydroDynamics* hydro; //!< Module for hydrodynamics.
	InputOutput* io; //!< Module for input/output.
	MPIHandler& mpihandler; //!< Handler for passing messages between processors.
	int steps; //!< Number of integration steps taken (used for switching order of operator split integration).

	//Structors.
	Integrator(IntegrationParameters& ipar, GridParameters& gpar, RadiationParameters& rpar,
			HydroParameters& hpar, PrintParameters& ppar, Scalings& scalings, MPIHandler& mpih);
	~Integrator();

	//Initialisation methods.
	void init();
	void init(std::string filename);

	//Integration methods.
	void march(const double& tmax);
	void march(const int& nsteps);

private:
	//Integration methods.
	double fluidStep();
	double calcTimeStep();
	double fluidStepSplitOld();
	double fluidStepSplit();
	void calcHydroFlux();
	void fluidStepSplitHydro(const double& dt);
	void fluidStepSplitRadiation(const double& dt);
	void advSolution(const double& dt);
	void fixSolution();
	void resetSrcTerms();
};



#endif /* INTEGRATOR_HPP_ */
