/** Provides the Hydrodynamics class.
 * @file Hydro.hpp
 *
 * @author Harrison Steggles
 * @date 13/01/2014 - the first version.
 * @date 16/01/2014 - modified code to accomodate new Boundary class.
 * @date 04/01/2014 - arguments now passed by const reference when appropriate.
 * @date 04/01/2014 - removed unnecessary variable declarations.
 * @date 21/02/2014 - globalUfromW and globalWfromU now uses memcpy rather than copying each element in a for loop.
 * @date 14/04/2014 - Hydro parameters now all held within a HydroParameters object.
 * @date 14/04/2014 - the function HLL has been added to solve the Riemann problem using the HLL scheme.
 * @date 14/04/2014 - new Riemann solver implemented. The function RotatedHLLC is now used to solve the Riemann
 * problem.
 * @date 01/05/2014 - calcFluxes now takes and integer argument: order, rather than using the fixed ORDER_S parameter in Grid3D.
 * @date 01/05/2014 - CFL now takes dt_max as argument instead of setting it from HydroParameters. It doesn't belong in HydroParameters.
 * @date 28/05/2014 - new methods added to update source terms instead of directly modifying the conservative fluid variables.
 * @date 21/07/2014 - replaced const reference to int/double with a const copy.
 * @date 24/11/2014 - now an Integrator subclass with lots of changes.
 */

#ifndef HYDRO_HPP_
#define HYDRO_HPP_

#include "Integrator.hpp"
#include "SlopeLimiter.hpp"
#include "Riemann.hpp"

class Fluid;
class HydroParameters;
class Radiation;
class Constants;

/**
 * @class HydroDynamics
 *
 * @brief Contains parameters and methods for hydrodynamic integration of a Fluid.
 *
 * @see Fluid
 * @see RiemannSolver
 * @see SlopeLimiter
 * @see Integrator
 *
 * @version 0.8, 24/11/2014
 */
class Hydrodynamics : public Integrator {

public:
	Hydrodynamics();
	~Hydrodynamics() {}

	void initialise(std::shared_ptr<Constants> c);

	virtual void preTimeStepCalculations(Fluid& fluid) const;
	virtual double calculateTimeStep(double dt_max, Fluid& fluid) const;
	virtual void integrate(double dt, Fluid& fluid) const;
	virtual void updateSourceTerms(double dt, Fluid& fluid) const;

	void setRiemannSolver(std::unique_ptr<RiemannSolver> riemannSolver);
	void setSlopeLimiter(std::unique_ptr<SlopeLimiter> slopeLimiter);

	//Calculation methods.
	void piecewiseLinear(FluidArray& Q_l, FluidArray& Q_c, FluidArray& Q_r, FluidArray& left_interp, FluidArray& right_interp) const;
	void reconstruct(Fluid& fluid) const;

	void calcFluxes(Fluid& fluid) const;
	//Integration methods.
	void updateBoundaries(Fluid& fluid) const;

private:
	std::shared_ptr<Constants> m_consts = nullptr;
	std::unique_ptr<RiemannSolver> m_riemannSolver = nullptr;
	std::unique_ptr<SlopeLimiter> m_slopeLimiter = nullptr;

	//Calculation methods.
	double soundSpeedSqrd(const double pre, const double den, const double gamma) const;
	double soundSpeed(const double pre, const double den, const double gamma) const;

	//Misc. methods.
	void Qisnan(const int id, const int i, const int xc, const int yc, const int zc) const;
};

#endif
