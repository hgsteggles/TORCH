/**
 * @file integrator.h
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
class HydroDynamics;
class Radiation;
class Thermodynamics;
class InputOutput;
class MPIHandler;
class IntegrationParameters;
class GridParameters;
class HydroParameters;
class RadiationParameters;
class PrintParameters;
class Scalings;

/**
 * @class Integrator
 *
 * @brief Class responsible for solving the state of a fluid.
 *
 * Integrator solves, in a step-wise manner, the state of a fluid system given parameters passed
 * into its constructor. The gas exists in a Grid3D object which is a collection of finite elements
 * called GridCell's. The hydrodynamics are solved using a rotated hybrid HLLC-HLL Riemann solver
 * to calculate fluxes on each GridCell face or GridJoin. Ionization from point source radiation
 * is implicitly solved and the column densities required for this are calculated via the method
 * of short characteristics. Heating/cooling from atomic processes is calculated using the
 * approximate functions in Henney (2009).
 *
 * @version 0.7, 13/06/2014
 */
class Integrator {
public:
	IntegrationParameters& iparams; //!< Contains integration parameters.
	Grid3D* grid; //!< The grid structure that holds integration cells.
	HydroDynamics* hydro; //!< Module for hydrodynamics.
	Radiation* rad; //!< Module for radiative transfer.
	Thermodynamics* thermo; //!< Module for thermodynamics.
	InputOutput* io; //!< Module for input/output.
	MPIHandler& mpihandler; //!< Handler for passing messages between processors.
	int steps; //!< Number of integration steps taken (used for switching order of operator split integration).

	//Structors.
	Integrator(IntegrationParameters& ipar, GridParameters& gpar, RadiationParameters& rpar,
			HydroParameters& hpar, PrintParameters& ppar, Scalings& scalings, MPIHandler& mpih);
	~Integrator();

	//Initialisation methods.
	void init(Scalings& scale);
	void init(std::string filename);

	//Integration methods.
	void march(const double& tmax);
	void march(const int& nsteps);

private:
	//Integration methods.
	double calcTimeStep();
	double fluidStepSplit();
	void updateBoundaries();
	void calcHydroFlux();
	void hydroFluidStepSplit(const double& dt);
	void radiationFluidStepSplit(const double& dt);
	void thermoFluidStepSplit(const double& dt);
	void advSolution(const double& dt);
	void fixSolution();
	void resetSrcTerms();
};



#endif /* INTEGRATOR_HPP_ */
