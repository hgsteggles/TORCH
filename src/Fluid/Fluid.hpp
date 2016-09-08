/**
 * Provides the Fluid class.
 * @file Fluid.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - The first version.
 * @date 03/12/2014 - Now correctly fixes solution in fixSolution().
 */

#ifndef FLUID_HPP_
#define FLUID_HPP_

#include <memory>

#include "Torch/Common.hpp"
#include "Torch/Parameters.hpp"
#include "Grid.hpp"
#include "Star.hpp"

class Constants;

/**
 * @class Fluid
 *
 * @brief The Fluid class contains a Grid and a Star, and provides methods for gas calculations.
 *
 * The Fluid should be passed into integrators so they have all the information needed to perform their integration. This class also provides methods for operating on every GridCell in the Grid.
 *
 * @see Grid
 * @see GridCell
 * @see Star
 */
class Fluid {
public:
	// Initialisers.
	void initialise(std::shared_ptr<Constants> c, FluidParameters fp);
	void initialiseGrid(GridParameters gp, StarParameters sp);

	// Updaters.
	void advSolution(const double dt);
	void fixSolution();
	void fixPrimitives();

	// Getters/Setters.
	Grid& getGrid();
	const Grid& getGrid() const;
	Star& getStar();
	const Star& getStar() const;

	// Conversion Methods.
	void globalWfromU();
	void globalUfromW();
	void globalQfromU();
	void globalUfromQ();

	// Calculations.
	double calcTemperature(double hii, double pre, double den) const;
	double calcSoundSpeed(double gamma, double pre, double den);
	double max(UID::ID id) const;
	double maxTemperature() const;
	double minTemperature() const;

	double heatCapacityRatio = 0;
	double massFractionH = 1.0; //!< Global mass fraction of hydrogen.
private:
	std::shared_ptr<Constants> consts = nullptr;
	Grid grid;
	Star star;
};

#endif // FLUID_HPP_
