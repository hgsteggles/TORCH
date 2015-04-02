#include "Boundary.hpp"

#include "Fluid/GridCell.hpp"
#include "MPI/MPI_Wrapper.hpp"

#include "GridCell.hpp"

/**
 * @brief Boundary constructor.
 * Requires a face and pointer to Grid3D so that the size of the 2D GridCell vector to be built is known.
 * @param face An integer specifying the face this Boundary is attached to (face%3 gives dimension, left
 * boundaries have face<3 and right boundaries have face>=3).
 * @param gptr Pointer to Grid3D instance that this Boundary is to be attached to.
 */

Boundary::Boundary(const int face, const int depth, Boundary::Type type)
	: m_face(face)
	, m_depth(depth)
	, m_type(type)
	, ncells(0)
	, m_lastGhostCellIterator(m_ghostCells.before_begin())
{

}

int Boundary::getFace() {
	return m_face;
}

Boundary::Type Boundary::getType() {
	return m_type;
}

int Boundary::getNumberGhosts() {
	return ncells;
}

std::forward_list<GridCell*> Boundary::getGhostCells() {
	return m_ghostCells;
}

GridCell* Boundary::addGhostCell(GridCell* copy_cell) {
	int dim = m_face%3;

	GridCell* bcptr = new GridCell();
	GridCell* newghosty = bcptr;
	for (int ig = 0; ig < m_depth - 1; ++ig) {
		GridCell* oldghost = newghosty;
		newghosty = new GridCell();

		if (m_face < 3)
			oldghost->left[dim] = newghosty;
		else
			oldghost->right[dim] = newghosty;
	}

	m_lastGhostCellIterator = m_ghostCells.insert_after(m_lastGhostCellIterator, bcptr);

	++ncells;
	return bcptr;
}

/**
 * @brief Boundary destructor.
 * Deletes all GridCell objects allocated
 */
Boundary::~Boundary() {
	int dim = m_face%3;
	for (GridCell* ghost : m_ghostCells) {
		GridCell* cptr = ghost;
		while (cptr != nullptr) {
			GridCell* oldcptr = cptr;
			if (m_face < 3)
				cptr = cptr->left[dim];
			else
				cptr = cptr->right[dim];
			delete oldcptr;
		}
	}
	m_ghostCells.clear();
}

/**
 * @brief ExternalBoundary constructor.
 * @param face
 * The face of the grid that the boundary is connected to (face < 3: left face, else: right face).
 * @param bcond
 * The boundary condition that must be applied to this boundary at the start of each hydro/radiation update step.
 * @param gptr
 * A pointer to the Grid3D object that the Boundary should be linked to with Grid3D::boundaryLink(Boundary*).
 */
ExternalBoundary::ExternalBoundary(const int face, const int depth, const Condition bcond)
	: Boundary(face, depth, Boundary::Type::EXTERNAL_BOUNDARY)
	, m_bc(bcond)
{

}

GridCell* ExternalBoundary::addGhostCell(GridCell* copy_cell) {
	GridCell* bcptr = Boundary::addGhostCell(copy_cell);

	int dim = m_face%3;

	bool requires_left_cell_wall = (m_face >= 3 && m_bc == Condition::PERIODIC) || (m_face < 3 && m_bc != Condition::PERIODIC);
	bool linking_left = (m_face >= 3);

	GridCell* ghostcell = bcptr;
	GridCell* icptr = copy_cell;

	while (ghostcell != nullptr) {
		if (linking_left) {
			ghostcell->left[dim] = icptr;
			ghostcell = ghostcell->right[dim];
		}
		else {
			ghostcell->right[dim] = icptr;
			ghostcell = ghostcell->left[dim];
		}

		if (m_bc == Condition::PERIODIC || m_bc == Condition::REFLECTING)
			icptr = requires_left_cell_wall ? icptr->right[dim] : icptr->left[dim];
	}

	return bcptr;
}

Condition ExternalBoundary::getCondition() {
	return m_bc;
}

/**
 * @brief Applies the boundary condition that was passed to this object through its constructor.
 */
void ExternalBoundary::applyBC() {
	int dim = m_face%3;
	for (GridCell* ghostcell : m_ghostCells) {
		for (GridCell *ghost = ghostcell, *cptr = nullptr; ghost != nullptr; ) {
			cptr  = (m_face < 3) ? ghost->right[dim] : ghost->left[dim];

			for(int iu = 0; iu < UID::N; iu++)
				ghost->Q[iu] = cptr->Q[iu];
			ghost->heatCapacityRatio = cptr->heatCapacityRatio;

			switch (m_bc) {
			case Condition::REFLECTING:
				ghost->Q[UID::VEL+dim] = -cptr->Q[UID::VEL+dim];
				break;
			case Condition::OUTFLOW:
				if((cptr->Q[UID::VEL+dim] > 0 && m_face < 3) || (cptr->Q[UID::VEL+dim] < 0 && m_face >= 3))
					ghost->Q[UID::VEL+dim] = -1.0*cptr->Q[UID::VEL+dim];
				break;
			case Condition::INFLOW:
				if((cptr->Q[UID::VEL+dim] < 0 && m_face < 3) || (cptr->Q[UID::VEL+dim] > 0 && m_face >= 3))
					ghost->Q[UID::VEL+dim] = -1.0*cptr->Q[UID::VEL+dim];
				break;
			case Condition::FREE:
				break;
			case Condition::PERIODIC:
				break;
			}

			ghost = (m_face < 3) ? ghost->left[dim] : ghost->right[dim];
		}
	}
	for (GridCell* ghostcell : m_ghostCells) {
		GridCell* ghost = ghostcell;
		GridCell* cptr = nullptr;
		if (m_face < 3)
			cptr = ghostcell->rjoin[dim]->rcell;
		else
			cptr = ghostcell->ljoin[dim]->lcell;
		for (int ir = 0; ir < RID::N; ++ir)
			ghost->R[ir] = cptr->R[ir];
	}
}

Partition::Partition(const int face, const int depth, const int dest)
	: Boundary(face, depth, Boundary::Type::PARTITION)
	, destination(dest)
	, nelements(1)
	, send_buffer(nullptr)
	, recv_buffer(nullptr)
{

}

Partition::~Partition() {
	if (send_buffer)
		delete[] send_buffer;
	if (recv_buffer)
		delete[] recv_buffer;
}

void Partition::load() {
	nelements = ncells*m_depth*(UID::N + 1);
	send_buffer = new double[nelements];
	recv_buffer = new double[nelements];
}

int Partition::getDestination() {
	return destination;
}

GridCell* Partition::addGhostCell(GridCell* copy_cell) {
	GridCell* bcptr = Boundary::addGhostCell(copy_cell);

	int dim = m_face%3;

	bool linking_left = (m_face >= 3);

	GridCell* ghostcell = bcptr;
	GridCell* icptr = copy_cell;

	while (ghostcell != nullptr) {
		if (linking_left) {
			ghostcell->left[dim] = icptr;
			ghostcell = ghostcell->right[dim];
		}
		else {
			ghostcell->right[dim] = icptr;
			ghostcell = ghostcell->left[dim];
		}

		icptr = linking_left ? icptr->left[dim] : icptr->right[dim];
	}

	return bcptr;
}

void Partition::applyBC() {
	MPIW& mpihandler = MPIW::Instance();
	int dim = m_face%3;
	int id = 0;

	for (GridCell* ghostcell : m_ghostCells) {
		for (GridCell *ghost = ghostcell, *cptr = nullptr; ghost != nullptr; ghost = (m_face < 3) ? ghost->left[dim] : ghost->right[dim]) {
			cptr  = (m_face < 3) ? ghost->right[dim] : ghost->left[dim];

			for(int iu = 0; iu < UID::N; ++iu)
				send_buffer[id++] = cptr->Q[iu];
			send_buffer[id++] = cptr->heatCapacityRatio;
		}
	}

	bool isPeriodic = (mpihandler.getRank() == 0 && m_face == 0) || (mpihandler.getRank() == mpihandler.nProcessors()-1 && m_face == 3);
	SendID tag = isPeriodic ? SendID::PERIODIC_MSG : SendID::PARTITION_MSG;

	mpihandler.exchange(send_buffer, recv_buffer, (int)nelements, destination, tag);

	id = 0;
	for (GridCell* ghostcell : m_ghostCells) {
		for (GridCell *ghost = ghostcell, *cptr = nullptr; ghost != nullptr; ghost = (m_face < 3) ? ghost->left[dim] : ghost->right[dim]) {
			for (int iu = 0; iu < UID::N; ++iu)
				ghost->Q[iu] = recv_buffer[id++];
			ghost->heatCapacityRatio = recv_buffer[id++];
		}
	}
}

void Partition::resetBuffer() {
	m_bufferCount = 0;
	m_recvCount = 0;
}

void Partition::addSendItem(double val) {
	if (m_bufferCount < (int)nelements)
		send_buffer[m_bufferCount++] = val;
}

double Partition::getRecvItem() {
	if (m_recvCount < m_bufferCount)
		return recv_buffer[m_recvCount++];
	else
		return 0;
}

void Partition::sendData(SendID tag) {
	MPIW::Instance().send(&m_bufferCount, 1, destination, tag);
	MPIW::Instance().send(send_buffer, m_bufferCount, destination, tag);
	resetBuffer();
}

void Partition::recvData(SendID tag) {
	resetBuffer();
	MPIW::Instance().receive(&m_bufferCount, 1, destination, tag);
	MPIW::Instance().receive(recv_buffer, m_bufferCount, destination, tag);
}
