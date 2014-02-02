/**
 * Provides the Star class.
 * @file star.hpp
 * @author Harrison Steggles
 * @date 29/01/2014, the first version.
 */

#ifndef STAR_HPP_
#define STAR_HPP_

/**
 * @class Star
 * @brief The star class holds information about a star's location and the nearest existing GridCell object to it on a
 * grid.
 * @version 0.3, 29/01/2014
 */
class Star {
public:
	GridCell* fcausal;
	int xc[3];
	int core;
	Star() : fcausal(NULL), core(0) {
		xc[0] = 0;
		xc[1] = 0;
		xc[2] = 0;
	}
	Star(int x, int y, int z) : fcausal(NULL), core(0) {
		xc[0] = x;
		xc[1] = y;
		xc[2] = z;
	}
	void setCore(int c) {
		core = c;
	}
};



#endif /* STAR_HPP_ */
