/**
 * Provides the Partition class (derived from the Boundary class).
 * @file partition.hpp
 * @author Harrison Steggles
 * @date 28/01/2014, the first version.
 */

#ifndef PARTITION_HPP_
#define PARTITION_HPP_

#include "boundary.hpp"
#include "mpihandler.hpp"
/**
 * @class Partition
 * @brief Boundary inherited class that acts as boundary between two grids handles by separate processors.
 * This class provides a Grid3D object with communication between processors so that processors can integrate
 * separate parts of the entire grid without erroneous results appearing at boundaries.
 * @version 0.3, 29/01/2014
 */
class Partition : public Boundary {
public:
	int destination;
	MPIHandler& mpihandler;
	Partition();
	Partition(int face, Grid3D* gptr, int dest, MPIHandler& mpih);
	void applyBC();
};



#endif /* PARTITION_HPP_ */
