/**
 * Provides the Boundary class.
 * @file boundary.hpp
 *
 * @author Harrison Steggles
 *
 * @date 16/01/2014 - the first version.
 * @date 04/01/2014 - arguments now passed by const reference when appropriate.
 * @date 08/05/2014 - integer specifying the number of ghost cell layers passed into Boundary constructor
 * now as it depends on spatial reconstruction information (ORDER_S) which now lies in IntegrationParameters
 * rather than GridParameters.
 */

#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include <vector>

class Grid3D;
class GridCell;

/**
 * @class Boundary
 * @brief Contains Grid3D face boundary cells and methods to apply boundary conditions during integration.
 * The Boundary class is abstract. Classes that inherit from Boundary will instantiate objects that contain GridCell objects that act as
 * "buffer" cells across a hydrodynamic grid face. Calling the Boundary::applyBC method will give any Boundary derived objects a chance to
 * set the buffer cells according to boundary conditions (ExternalBoundary) or to receive information from processors processing adjacent
 * grids (Partition).
 * Grid3D handles the binding of a Boundary to its grid of GridCell objects.
 * @see ExternalBoundary
 * @see Partition
 * @see Grid3D
 * @see GridCell
 *
 * @version 0.4, 04/02/2014
 */
class Boundary {
public:
	std::vector<std::vector<GridCell*> > ghostcells; //!< 2D vector of buffer cells that cover the face of a grid.
	int face; //!< An id for the grid face this Boundaryis attached to.
	int nghosts; //!< Number of buffer layers. Depends on temporal and spatial order of integration.
	bool isPartition; //!< Is this Boundary instance a Partition?

	Boundary(const int& face, const int nOfGhosts, Grid3D* gptr);
	virtual ~Boundary();

	virtual void applyBC() = 0;
};



#endif /* BOUNDARY_HPP_ */
