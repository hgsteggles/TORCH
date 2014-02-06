/**
 * @file integrator.hpp
 *
 * @author Harrison Steggles
 * @date 04/02/2014 - initial verision. Encapsulating program into integrator class.
 */

#ifndef INTEGRATOR_HPP_
#define INTEGRATOR_HPP_

#include <string>

class Grid3D;
class Radiation;
class HydroDynamics;
class InputOutput;
class MPIHandler;
class GridParameters;
class RadiationParameters;
class HydroParameters;
class PrintParameters;
class Scalings;

class Integrator {
public:
	Integrator(GridParameters& gpar, RadiationParameters& rpar,
			HydroParameters& hpar, PrintParameters& ppar, Scalings& scalings, MPIHandler& mpih);
	~Integrator();
	void init();
	void init(std::string filename);
	/**
	 * @brief Solves the state of radiation and hydrodynamic variables until the time tmax.
	 * Marches the conserved variables, Grid3D::U, in a Grid3D object until unitless (simulation) time tmax has been reached.
	 * @param tmax The time (unitless) that march stops marching.
	 */
	void march(const double& tmax);
	/**
	 * @brief Solves the state of radiation and hydrodynamic variables nsteps times.
	 * Marches the conserved variables, Grid3D::U, in a Grid3D object nsteps times.
	 * @param nsteps The number of steps to march.
	 */
	void march(const int& nsteps);
	double fluidStep();

	Grid3D* grid;
	Radiation* rad;
	HydroDynamics* hydro;
	InputOutput* io;
	MPIHandler& mpihandler;
};



#endif /* INTEGRATOR_HPP_ */
