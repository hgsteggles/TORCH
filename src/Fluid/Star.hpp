/**
 * Provides the Star class.
 *
 * @file Star.hpp
 *
 * @author Harrison Steggles
 *
 * @date 29/01/2014 - Version 0.3 pushed to Master.
 * @date 30/01/2014 - Changed type of Star::xc to double, location now in units of cell width and origin at corner
 * of grid NOT at centre of GridCell in corner.
 * @date 30/01/2014 - Changed type of Star::xc back to int but coordinate system stays the same so that Star position
 * can snap to vertex locations.
 * @date 21/07/2014 - Constructor now receives a StarParameters and Grid3D reference so it can build itself (the
 * Radiation class was doing this before). Also, Star holds a reference to StarParameters instead of copying StarParameters'
 * members.
 * @date 24/11/2014 - added wind parameters and injection methods.
 */

#ifndef STAR_HPP_
#define STAR_HPP_

#include <array>
#include <memory>

#include "Torch/Common.hpp"
#include "Torch/Parameters.hpp"

class GridCell;
class Grid;
class Constants;

/**
 * @class Star
 *
 * @brief The star class contains stellar properties and methods for injecting a stellar wind.
 *
 * @see StarParameters
 * @see GridCell
 * @see Grid
 */
class Star {
public:
	enum class Location : unsigned int { LEFT, HERE, RIGHT };
	std::array<double, 3> xc = std::array<double, 3>{{ 0, 0, 0 }};

	void initialise(std::shared_ptr<Constants> c, StarParameters sp, Location containing_core, const Vec3& delta_x);

	void setWindCells(Grid& grid);

	void injectEnergyMomentum(Grid& grid);

	void fixDensityPressure(Grid& grid);

	int getWindCellRadius() const;

	bool on = false;
	Vec3 dx =  Vec3{{ 0, 0, 0 }};
	double photonEnergy = 0;
	double photonRate = 0;
	double massLossRate = 0;
	double windVelocity = 0;
	double windTemperature = 0;
	int windCellRadius = 0;
	Location core = Location::HERE;

private:
	std::shared_ptr<Constants> consts = nullptr;
	double mdot = 0;
	double edot = 0;
};

#endif // STAR_HPP_
