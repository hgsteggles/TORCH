/**
 * Provides the Boundary class.
 * @file boundary.hpp
 *
 * @author Harrison Steggles
 *
 * @date 16/01/2014, the first version.
 */

#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include <vector>
#include "parameters.hpp"
#include "grid3d.hpp"

class Grid3D;
/**
 * @class Boundary
 * @brief Holds a Grid3D face boundary to apply boundary conditions during integration.
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
 * @version 0.3, 29/01/2014
 */
class Boundary {
public:
	std::vector<std::vector<GridCell*> > ghostcells; //!< 2D vector of buffer cells that cover the face of a grid.
	int face; //!< An id for the grid face this Boundaryis attached to.
	int nghosts; //!< Number of buffer layers. Depends on temporal and spatial order of integration.
	bool isPartition; //!< Is this Boundary instance a Partition?
	/**
	 * @brief Boundary constructor.
	 * Requires a face and pointer to Grid3D so that the size of the 2D GridCell vector to be built is known.
	 * @param face The face this Boundary is attached to.
	 * @param gptr Pointer to Grid3D instance that this Boundary is to be attached to.
	 */
	Boundary(int face, Grid3D* gptr);
	/**
	 * @brief Boundary destructor.
	 */
	virtual ~Boundary();
	/**
	 * @brief Applies changes to Boundary cells.
	 * Boundary cells are modified according to boundary conditions or information from other processors.
	 */
	virtual void applyBC() = 0;
};



#endif /* BOUNDARY_HPP_ */
