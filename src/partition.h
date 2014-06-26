/**
 * Provides the Partition class (derived from the Boundary class).
 * @file partition.h
 * @author Harrison Steggles
 * @date 28/01/2014, the first version.
 * @date 04/01/2014, arguments now passed by const reference when appropriate.
 * @date 21/02/2014 - uses new MPIHandler function called exchange which exchanges data with neighbouring processes.
 * A copy of the send data is made when receiving before sending (blocking communication on odd processors).
 */

#ifndef PARTITION_HPP_
#define PARTITION_HPP_

#include "boundary.h"

class MPIHandler;
class Grid3D;

/**
 * @class Partition
 *
 * @brief Boundary inherited class that acts as boundary between two grids handles by separate processors.
 *
 * This class provides a Grid3D object with communication between processors so that processors can integrate
 * separate parts of the entire grid without erroneous results appearing at boundaries.
 *
 * @version 0.7, 13/06/2014
 */
class Partition : public Boundary {
public:
	int destination;
	MPIHandler& mpihandler;
	Partition();
	Partition(const int face, const int nOfGhosts, Grid3D* gptr, const int& dest, MPIHandler& mpih);
	void applyBC();
};



#endif /* PARTITION_HPP_ */
