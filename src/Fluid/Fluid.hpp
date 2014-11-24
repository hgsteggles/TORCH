/**
 * Provides the Fluid class.
 * @file Fluid.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef FLUID_HPP_
#define FLUID_HPP_

#include "Parameters.hpp"
#include "Star.hpp"
#include "Grid.hpp"
#include "GridCell.hpp"

#include <memory>
#include <cmath>

class Constants;

/**
 * @class Fluid
 *
 * @brief The Fluid class contains a Grid and a Star and provides methods for gas calculations.
 *
 * The Fluid should be passed into integrators so they have all the information needed to perform their integration. This class
 * also provides methods for operating on every GridCell in the Grid.
 *
 * @see Grid
 * @see GridFactory
 * @see GridCell
 * @see Star
 *
 * @version 0.8, 24/11/2014
 */
class Fluid {
public:
	void initialise(std::shared_ptr<Constants> c, FluidParameters fp);

	void initialiseGrid(GridParameters gp, StarParameters sp);

	// Conversion Methods.
	void globalWfromU() const;
	void globalUfromW() const;
	void globalQfromU() const;
	void globalUfromQ() const;

	void advSolution(const double dt);
	void fixSolution();
	void fixPrimitives();

	// Calculations.
	double calcTemperature(double hii, double pre, double den);
	double calcSoundSpeed(double gamma, double pre, double den);

	Grid& getGrid();
	const Grid& getGrid() const;
	Star& getStar();
	const Star& getStar() const;

	double heatCapacityRatio = 0; //!<
	double massFractionH = 1.0;
private:
	std::shared_ptr<Constants> consts = nullptr;
	Grid grid;
	Star star;
};

#endif /* FLUID_HPP_ */
