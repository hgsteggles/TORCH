/** Provides the Integrator abstract base class.
 * @file Integrator.hpp
 *
 * @author Harrison Steggles
 * @date 24/11/2014 - the first version.
 */

#ifndef INTEGRATOR_HPP_
#define INTEGRATOR_HPP_

#include <string>

class Fluid;

/**
 * @class Integrator
 *
 * @brief Abstract base class for integration of Fluid for a physics sub-problem.
 *
 * @see Fluid
 */
class Integrator {
public:
	Integrator(std::string name) : componentName(name) { }
	virtual ~Integrator() { }
	virtual void preTimeStepCalculations(Fluid& fluid) const = 0;
	virtual double calculateTimeStep(double dt_max, Fluid& fluid) const = 0;
	virtual void integrate(double dt, Fluid& fluid) const = 0;
	virtual void updateSourceTerms(double dt, Fluid& fluid) const = 0;

	virtual std::string getComponentName() final {
		return componentName;
	}
	std::string componentName = "DefaultComponentName";
};



#endif // INTEGRATOR_HPP_
