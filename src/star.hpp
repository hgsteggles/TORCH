/**
 * Provides the Star class.
 * @file star.hpp
 * @author Harrison Steggles
 * @date 29/01/2014 - Version 0.3 pushed to Master.
 * @date 30/01/2014 - Changed type of Star::xc to double, location now in units of cell width and origin at corner of grid NOT at centre of GridCell in corner.
 * @date 30/01/2014 - Changed type of Star::xc back to int but coordinate system stays the same so that Star position can snap to vertex locations.
 */

#ifndef STAR_HPP_
#define STAR_HPP_

#include <stddef.h>

class GridCell;

/**
 * @class Star
 * @brief The star class holds information about a star's location and the nearest existing GridCell object to it on a grid.
 * @version 0.4, 30/01/2014
 */
class Star {
public:
	GridCell* fcausal;
	int x[3];
	int core;
	Star() : fcausal(NULL), core(0) {
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
	}
	Star(int x, int y, int z) : fcausal(NULL), core(0) {
		this->x[0] = x;
		this->x[1] = y;
		this->x[2] = z;
	}
	void setCore(int c) {
		core = c;
	}
};



#endif /* STAR_HPP_ */
