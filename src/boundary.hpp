/*
 * boundary.hpp
 *
 *  Created on: 14 Jan 2014
 *      Author: harry
 */

#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include <vector>
#include "parameters.hpp"
#include "grid3d.hpp"
class Grid3D;

class Boundary {
public:
	Boundary(int face, Condition bc, Grid3D* gptr);
	~Boundary();
	void applyBC();
	std::vector<std::vector<GridCell*> > ghostcells;
	int face, nghosts;
	Condition bc;
};



#endif /* BOUNDARY_HPP_ */
