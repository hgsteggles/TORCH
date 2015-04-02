/**
 * Provides the Boundary class and derived classes.
 * @file Boundary.hpp
 *
 * @author Harrison Steggles
 *
 * @date 16/01/2014 - The first version.
 * @date 04/01/2014 - Arguments now passed by const reference when appropriate.
 * @date 08/05/2014 - Integer specifying the number of ghost cell layers passed into Boundary constructor
 * now as it depends on spatial reconstruction information (ORDER_S) which now lies in IntegrationParameters
 * rather than GridParameters.
 * @date 21/07/2014 - Replaced const reference to int/double with a const copy.
 * @date 24/11/2014 - Now add GridCells to Boundary by calling Boundary::addGhostCell. Passed in GridCell gets linked so that BCs can
 * be quickly calculated.
 * @date 24/11/2014 - Derived classes: Partition and ExternalBoundary have been moved to here.
 * @date 24/11/2014 - Partition now holds send AND receive buffers for sending data across processors rather than allocating heap
 * space for a buffer at every call to MPIW::exchange (which has also been modified).
 * @date 15/12/2014 - Fixed bug in ExternalBoundary::applyBC where heatCapacityRatio not sent to ghost cells.
 */

#ifndef BOUNDARY_HPP_
#define BOUNDARY_HPP_

#include <forward_list>

#include "MPI/MPI_Wrapper.hpp"
#include "Torch/Common.hpp"

class GridCell;

/**
 * @class Boundary
 *
 * @brief Contains Grid face boundary cells and methods to apply boundary conditions during integration.
 *
 * The Boundary class is abstract. Classes that inherit from Boundary will instantiate objects that contain GridCell objects that act as
 * "buffer" cells across a hydrodynamic grid face. Calling the Boundary::applyBC method will give any Boundary derived objects a chance to
 * set the buffer cells according to boundary conditions (ExternalBoundary) or to receive information from processors processing adjacent
 * Grids (Partition).
 * GridFactory handles the binding of a Boundary to a Simulation Grid of GridCells.
 *
 * @see ExternalBoundary
 * @see Partition
 * @see Grid
 * @see GridFactory
 * @see GridCell
 *
 * @version 0.8, 24/11/2014
 */
class Boundary {
public:
	enum class Type : unsigned int {PARTITION, EXTERNAL_BOUNDARY, NONE};

	Boundary(const int face, const int depth, Boundary::Type type);
	virtual ~Boundary();

	int getFace();
	Boundary::Type getType();
	int getNumberGhosts();
	virtual GridCell* addGhostCell(GridCell* copy_cell);
	std::forward_list<GridCell*> getGhostCells();

	virtual void applyBC() = 0;

protected:
	int m_face; //!< An id for the grid face this Boundaryis attached to.
	int m_depth; //!< Number of buffer layers. Depends on temporal and spatial order of integration.
	Boundary::Type m_type;
	int ncells;
	std::forward_list<GridCell*> m_ghostCells;
	std::forward_list<GridCell*>::iterator m_lastGhostCellIterator;
};

/**
 * @class ExternalBoundary
 *
 * @brief Holds a Grid face ExternalBoundary to apply boundary conditions during integration.
 *
 * The ExternalBoundary class inherits from the Boundary class. ExternalBoundary class contain GridCells that buffer
 * a Grid face from unsimulated regions. ExternalBoundary::applyBC applies boundary conditions to the contained GridCells.
 * GridFactory handles the binding of an ExternalBoundary to a grid of GridCells.
 *
 * @see Boundary
 * @see Grid
 * @see GridFactory
 * @see GridCell
 *
 * @version 0.8, 24/11/2014
 */
class ExternalBoundary : public Boundary {
public:
	ExternalBoundary(const int face, const int depth, const Condition bcond);

	Condition getCondition();
	GridCell* addGhostCell(GridCell* copy_cell);
	void applyBC();

private:
	Condition m_bc;
};

/**
 * @class Partition
 *
 * @brief Boundary inherited class that acts as boundary between two Grids handled by separate processors.
 *
 * This class provides a Grid with communication between processors so that processors can integrate
 * separate parts of the entire Grid without erroneous results appearing at boundaries.
 *
 * @see Boundary
 * @see Grid
 * @see GridFactory
 * @see GridCell
 *
 * @version 0.8, 24/11/2014
 */
class Partition : public Boundary {
public:
	Partition();
	Partition(const int face, const int depth, const int dest);
	~Partition();

	void load();
	GridCell* addGhostCell(GridCell* copy_cell);

	int getDestination();
	void applyBC();

	void resetBuffer();
	void addSendItem(double val);
	double getRecvItem();
	void sendData(SendID tag);
	void recvData(SendID tag);
private:
	int destination = 0;
	unsigned int nelements = 0;
	int m_bufferCount = 0;
	int m_recvCount = 0;
	int m_sendCount = 0;
	double* send_buffer = nullptr;
	double* recv_buffer = nullptr;
};

#endif // BOUNDARY_HPP_
