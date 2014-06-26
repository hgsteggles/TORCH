/**
 * Provides the ExternalBoundary class.
 * @file externalboundary.h
 *
 * @author Harrison Steggles
 *
 * @date 29/01/2014 - the first version.
 * @date 04/01/2014 - arguments now passed by const reference when appropriate.
 * @date 05/02/2014 - solved buffer overflow problem for non-reflecting boundary conditions - replaced 'j' variable
 * in Q array with 'iu'.
 * @date 21/02/2014 - applyBC now copies grid face radiation properties (cell and path optical depths) onto first
 * layer of ghost cells for
 * raytracing from a star that has snapped to vertex that lies on edge of grid.
 * @date 08/05/2014 - integer specifying the number of ghost cell layers passed into ExternalBoundary constructor
 * now as it depends on spatial reconstruction information (ORDER_S) which now lies in IntegrationParameters
 * rather than GridParameters.
 */

#ifndef EXTERNAL_HPP_
#define EXTERNAL_HPP_

#include "boundary.h"
#include "constants.h"

class Grid3D;

/**
 * @class ExternalBoundary
 *
 * @brief Holds a Grid3D face ExternalBoundary to apply boundary conditions during integration.
 *
 * The ExternalBoundary class inherits from the Boundary class. ExternalBoundary instances provide GridCell instances
 * that buffer a grid face from unsimulated regions. ExternalBoundary::applyBC applies boundary conditions to the
 * contained GridCell instances. Grid3D handles the binding of an ExternalBoundary to its grid of GridCell objects.
 *
 * @see Boundary
 * @see Grid3D
 * @see GridCell
 *
 * @version 0.7, 13/06/2014
 */
class ExternalBoundary : public Boundary {
public:
	Condition bc;

	ExternalBoundary(const int face, const int nOfGhosts, const Condition& bcond, Grid3D* gptr);

	void applyBC();
};



#endif /* EXTERNAL_HPP_ */
