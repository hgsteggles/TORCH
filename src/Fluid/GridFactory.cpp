#include "GridFactory.hpp"
#include "GridCell.hpp"
#include "Boundary.hpp"
#include "Parameters.hpp"
#include "ProgressBar.hpp"

#include <iostream>
#include <cmath>

static int calcCoreCells(int nx, int nproc, int iproc) {
	int ncells = nx/nproc;
	int cells_left = nx - nproc*ncells;
	if (nproc - iproc <= cells_left)
		++ncells;
	return ncells;
}

static int calcLeftBoundaryPosition(int nx, int nproc, int iproc) {
	int start_xc = 0;
	for (int oproc = 0; oproc < iproc; ++oproc)
		start_xc += calcCoreCells(nx, nproc, oproc);
	return start_xc;
}

Grid GridFactory::createGrid(const std::shared_ptr<Constants>& c, const GridParameters& gp, const std::array<int, 3>& star_pos, int radius) {
	MPIW& mpihandler = MPIW::Instance();
	int rproc = mpihandler.getRank();
	int nproc = mpihandler.nProcessors();

	// Set up Grid parameters.
	Grid grid;
	grid.ncells = gp.ncells;
	grid.sideLength = gp.sideLength;
	grid.spatialOrder = gp.spatialOrder;
	grid.geometry = c->geometryParser.parseEnum(gp.geometry);
	for (int i = 0; i < c->nd; ++i)
		grid.dx[i] = gp.sideLength/(double)gp.ncells[0];

	grid.coreCells = gp.ncells;
	grid.coreCells[0] = calcCoreCells(gp.ncells[0], mpihandler.nProcessors(), mpihandler.getRank());
	if (grid.coreCells[0] == 0) {
		throw std::runtime_error("GridBuilder::createGrid: zero cells in processor (" +  std::to_string(mpihandler.getRank()) + ").");
		return grid;
	}

	// Build the doubly linked list.
	int start_xc = calcLeftBoundaryPosition(grid.ncells[0], nproc, rproc);
	GridCellPair firstLastCell = buildGrid(grid.coreCells, c->nd, start_xc, grid.m_joins);

	// Zip the GridCells into a DualIntrusiveContainer.
	grid.m_cells = std::move(zipList(firstLastCell, star_pos, radius, c->nd));

	// Calculate volumes.
	for (GridCell& cell : grid.m_cells.getAll())
		cell.vol = vol_cell(cell.xc[0], grid.dx, grid.geometry, c->nd);

	// Parse boundary condition into enum.
	ConditionVecPair leftRightBC;
	for (unsigned int i = 0; i < 3; ++i) {
		leftRightBC.first[i] = c->conditionParser.parseEnum(gp.leftBC[i]);
		leftRightBC.second[i] = c->conditionParser.parseEnum(gp.rightBC[i]);
	}

	// Build the boundaries.
	BoundaryVecPair boundaryPair = buildBoundaries(leftRightBC, gp.spatialOrder, c->nd);
	grid.m_leftBoundaries = boundaryPair.first;
	grid.m_rightBoundaries = boundaryPair.second;

	// Link the boundaries to the Grid.
	for (int i = 0; i < 3; ++i) {
		boundaryLink(grid.m_leftBoundaries[i], firstLastCell.first, grid.coreCells, grid.m_joins);
		boundaryLink(grid.m_rightBoundaries[i], firstLastCell.first, grid.coreCells, grid.m_joins);
	}

	// Calculate Join areas.
	for (int i = 0; i < c->nd; ++i) {
		for (GridJoin& join : grid.m_joins[i])
			join.area = area_join(join.xj, i, grid.dx, grid.geometry, c->nd);
	}

	// Set Grid wall locations.
	grid.setLeftX((int)firstLastCell.first->xc[0]);
	grid.setRightX((int)firstLastCell.second->xc[0]);

	return std::move(grid);
}

GridFactory::GridCellPair GridFactory::buildGrid(const Coords& ncells, int nd, int start_xc, std::array<JoinContainer, 3>& joins) {
	int proc_rank = MPIW::Instance().getRank();
	GridCell* firstCell = nullptr;
	GridCell* lastCell = nullptr;

	ProgressBar progBar = ProgressBar(100, 5, "Building Grid3D", false);

	for (int i = 0; i < ncells[0]*ncells[1]*ncells[2]; ++i){
		double dummy;
		progBar.update((100*i/(double)(ncells[0]*ncells[1]*ncells[2])), dummy, proc_rank == 0);
		GridCell* newcptr = new GridCell();
		// Set GridCell index and link to left and right neighbours.
		if (lastCell != nullptr) {
			if ((int)lastCell->xc[0] != (int)firstCell->xc[0] + ncells[0] - 1) {
				newcptr->set_xcs(lastCell->xc[0] + 1, lastCell->xc[1], lastCell->xc[2]);
				joins[0].push_back(link(0, lastCell, newcptr));
				joins[1].push_back(link(1, Grid::locate(Coords{(int)newcptr->xc[0], (int)newcptr->xc[1] - 1, (int)newcptr->xc[2]}, firstCell), newcptr));
				joins[2].push_back(link(2, Grid::locate(Coords{(int)newcptr->xc[0], (int)newcptr->xc[1], (int)newcptr->xc[2] - 1}, firstCell), newcptr));
			}
			else if ((int)lastCell->xc[0] == (int)firstCell->xc[0] + ncells[0] - 1 && (int)lastCell->xc[1] != ncells[1] - 1) {
				newcptr->set_xcs(firstCell->xc[0], lastCell->xc[1] + 1, lastCell->xc[2]);
				joins[1].push_back(link(1, Grid::locate(Coords{(int)newcptr->xc[0], (int)newcptr->xc[1] - 1, (int)newcptr->xc[2]}, firstCell), newcptr));
				joins[2].push_back(link(2, Grid::locate(Coords{(int)newcptr->xc[0], (int)newcptr->xc[1], (int)newcptr->xc[2] - 1}, firstCell), newcptr));
			}
			else if ((int)lastCell->xc[0] == (int)firstCell->xc[0] + ncells[0] - 1 && (int)lastCell->xc[1] == ncells[1] - 1) {
				newcptr->set_xcs(firstCell->xc[0], firstCell->xc[1], lastCell->xc[2] + 1);
				joins[2].push_back(link(2, Grid::locate(Coords{(int)newcptr->xc[0], (int)newcptr->xc[1], (int)newcptr->xc[2]-1}, firstCell), newcptr));
			}
		}
		else {
			newcptr->xc[0] = start_xc + 0.5;
			for (unsigned int i = 1; i < nd; ++i)
				newcptr->xc[i] = 0.5;
			firstCell = newcptr;
		}
		lastCell = newcptr;
	}

	progBar.end(proc_rank == 0);
	return std::make_pair(firstCell, lastCell);
}

DualIntrusiveContainer<GridCell> GridFactory::zipList(const GridCellPair& firstLastCell, const Coords& xs, int radius, int nd) {
	GridCell* fcausal = Grid::getNearestCell(xs, firstLastCell);
	GridCell *oldwptr = nullptr, *startcptr = nullptr, *oldcptr = nullptr;
	GridCell *startwptr = Grid::locate(xs, firstLastCell.first);

	CellContainer windCells, causalCells;

	for (GridCell *cptr = fcausal; cptr != nullptr; cptr = Grid::nextCausal(cptr, fcausal, nd)) {
		double dist2 = 0;
		for (unsigned int idim = 0; idim < nd; ++idim)
			dist2 += (cptr->xc[idim] - xs[idim])*(cptr->xc[idim] - xs[idim]);
		if (dist2 <= radius*radius || cptr == startwptr)
			windCells.push_back(cptr);
		else
			causalCells.push_back(cptr);
	}

	return DualIntrusiveContainer<GridCell>(windCells, causalCells);
}

std::pair<BoundaryVec, BoundaryVec> GridFactory::buildBoundaries(const ConditionVecPair& leftRightBC, int spatialOrder, int nd) {
	MPIW& mpihandler = MPIW::Instance();
	int nproc = mpihandler.nProcessors();
	int rproc = mpihandler.getRank();

	BoundaryVec leftBoundaries, rightBoundaries;
	const ConditionVec& leftBC = leftRightBC.first;
	const ConditionVec& rightBC = leftRightBC.second;

	int left_rank = (rproc - 1 < 0) ? rproc - 1 + nproc : rproc - 1;
	int right_rank = (rproc + 1 >= nproc) ? rproc + 1 - nproc : rproc + 1;
	bool periodic_and_mpi = nproc > 1 && leftBC[0] == Condition::PERIODIC && rightBC[0] == Condition::PERIODIC;

	for (int i = 0; i < 3; ++i)
		leftBoundaries[i] = rightBoundaries[i] = nullptr;

	if (nd > 0) {
		if (mpihandler.getRank() != 0 || periodic_and_mpi)
			leftBoundaries[0] = new Partition(0, spatialOrder+1, left_rank);
		else
			leftBoundaries[0] = new ExternalBoundary(0, spatialOrder+1, leftBC[0]);
		if (rproc != nproc - 1 || periodic_and_mpi)
			rightBoundaries[0] = new Partition(3, spatialOrder+1, right_rank);
		else
			rightBoundaries[0] = new ExternalBoundary(3, spatialOrder+1, rightBC[0]);
	}
	if (nd > 1) {
		leftBoundaries[1] = new ExternalBoundary(1, spatialOrder+1, leftBC[1]);
		rightBoundaries[1] = new ExternalBoundary(4, spatialOrder+1, rightBC[1]);
	}
	if (nd > 2) {
		leftBoundaries[2] = new ExternalBoundary(2, spatialOrder+1, leftBC[2]);
		rightBoundaries[2] = new ExternalBoundary(5, spatialOrder+1, rightBC[2]);
	}

	return std::make_pair(leftBoundaries, rightBoundaries);
}

/**
 * @brief Double-links the passed GridCell objects to a GridJoin object.
 * Links two GridCell objects to a common GridJoin object via GridCell::rjoin and GridCell::ljoin.
 * GridJoin contains links to both GridCell objects via GridJoin::rcell and GridJoin::lcell. The
 * coordinates of the GridJoin object is set according to those of the GridCell objects.
 * @param dim
 * Dimension across which to set up links.
 * @param lcptr
 * Pointer to left GridCell object.
 * @param rcptr
 * Pointer to right GridCell object.
 */
GridJoin* GridFactory::weakLink(const int dim, GridCell* lcptr, GridCell* rcptr){
	if (lcptr != nullptr && rcptr != nullptr) {
		GridJoin* jptr = new GridJoin();
		rcptr->ljoin[dim] = jptr;
		lcptr->rjoin[dim] = jptr;
		jptr->rcell = rcptr;
		jptr->lcell = lcptr;
		for(int i = 0; i < 3; ++i)
			jptr->xj[i] = jptr->rcell->xc[i];
		jptr->xj[dim] = jptr->rcell->xc[dim] - 0.5;
		return jptr;
	}
	else
		return nullptr;
}

/**
 * @brief Double-links the passed GridCell objects to each other and a GridJoin object.
 * Links together two GridCell objects via GridCell::left and GridCell::right and also links them
 * both to a common GridJoin object via GridCell::rjoin and GridCell::ljoin. GridJoin contains links
 * to both GridCell objects via GridJoin::rcell and GridJoin::lcell. The coordinates of the
 * GridJoin object is set according to those of the GridCell objects.
 *
 * @param dim
 * Dimension across which to set up links.
 * @param lcptr
 * Pointer to left GridCell object.
 * @param rcptr
 * Pointer to right GridCell object.
 */
GridJoin* GridFactory::link(const int dim, GridCell* lcptr, GridCell* rcptr){
	if (lcptr != nullptr && rcptr != nullptr) {
		rcptr->left[dim] = lcptr;
		lcptr->right[dim] = rcptr;
		GridJoin* jptr = weakLink(dim, lcptr, rcptr);
		return jptr;
	}
	else
		return nullptr;
}

/**
 * @brief Links GridCell objects on a grid face with GridCell objects in a Boundary object via a
 * GridJoin object. Each GridCell object on a grid face is linked via a GridJoin to a GridCell object
 * owned by a Boundary. The GridCell objects are not double-linked to each other because traversal
 * over the simulation grid is achieved by detecting nullptr pointers on the boundary.
 * @param bptr
 * Pointer to the Boundary object to be linked.
 */
void GridFactory::boundaryLink(Boundary* bptr, GridCell* const firstCell, const Coords& coreCells, std::array<JoinContainer, 3>& joins) {
	if (bptr != nullptr) {
		int dim = bptr->getFace()%3;
		bool isLeft = bptr->getFace() < 3;

		GridCell* link_cell = isLeft ? firstCell : Grid::traverse1D(dim, coreCells[dim]-1, firstCell);
		GridCell* copy_cell = nullptr;

		bool isPeriodic = false;
		if (bptr->getType() == Boundary::Type::EXTERNAL_BOUNDARY) {
			ExternalBoundary* ebptr = static_cast<ExternalBoundary*>(bptr);
			isPeriodic = ebptr->getCondition() == Condition::PERIODIC;
		}
		bool req_left_face = (!isLeft && isPeriodic) || (isLeft && !isPeriodic);
		copy_cell = req_left_face ? firstCell : Grid::traverse1D(dim, coreCells[dim]-1, firstCell);

		while (link_cell != nullptr) {
			GridCell* bcptr = bptr->addGhostCell(copy_cell);

			for (int id = 0; id < 3; ++id) {
				bcptr->xc[id] = link_cell->xc[id];
			}
			bcptr->xc[dim] += isLeft ? -1 : 1;


			joins[dim].push_back( weakLink(dim, isLeft ? bcptr : link_cell, isLeft ? link_cell : bcptr) );


			link_cell = Grid::nextCell2D(dim, link_cell);
			if (copy_cell != nullptr)
				copy_cell = Grid::nextCell2D(dim, copy_cell);
		}

		if (bptr->getType() == Boundary::Type::PARTITION) {
			Partition* partition = static_cast<Partition*>(bptr);
			partition->load();
		}
	}
}

/**
 * @brief Calculates the area between two GridCell objects given its location and the geometry of Grid3D.
 * @param xj
 * Coordinates of the GridJoin between the two GridCell objects.
 * @param dim
 * Direction normal to the face of the GridJoin.
 * @return The area.
 */
double GridFactory::area_join(const Vec3& xj, const int dim, const Vec3& dx, Geometry geometry, int nd) {
	double PI = 3.14159265359;
	double area = 0, r1, r2;
	if (geometry == Geometry::CARTESIAN) {
		area = 1.0;
		for (int i = 0; i < nd; ++i) {
			if (i != dim)
				area *= dx[i];
		}
	}
	else if (geometry == Geometry::CYLINDRICAL){
		/* Cylindrical polars (ND <= 2) */
		if (dim == 0) {
			r1 = xj[0]*dx[0];
			area = 2.0*PI*r1;
			if (nd >= 2)
				area *= dx[1];
		}
		else if (dim == 1) {
			r1 = xj[0]*dx[0] - 0.5*dx[0];
			r2 = xj[0]*dx[0] + 0.5*dx[0];
			area = PI*(r2 - r1)*(r2 + r1);
		}
	}
	else if (geometry == Geometry::SPHERICAL) {
		/* Spherical polars (ND = 1 only) */
		r1 = xj[0]*dx[0];
		area = 4.0*PI*r1*r1;
	}
	return area;
}

/**
 * @brief Calculates the volume of a GridCell given its location and the geometry of Grid3D.
 * @param cptr
 * The GridCell object for which the volume will be calculated.
 */
double GridFactory::vol_cell(double rc, const Vec3& dx, Geometry geometry, int nd) {
	double PI = 3.14159265359;
	double volume = 0, r1, r2;
	if (geometry == Geometry::CARTESIAN) {
		/* Cartesian */
		volume = 1.0;
		for (int i = 0; i < nd; ++i)
			volume *= dx[i];
	}
	else if (geometry == Geometry::CYLINDRICAL) {
		/* Cylindrical polars (ND <= 2) */
		r2 = rc*dx[0] + 0.5*dx[0];
		r1 = rc*dx[0] - 0.5*dx[0];
		volume = PI*(r2 - r1)*(r2 + r1);
		if (nd >= 2)
			volume *= dx[1];
	}
	else if (geometry == Geometry::SPHERICAL) {
		/* Spherical polars (ND = 1 only) */
		r2 = rc*dx[0] + 0.5*dx[0];
		r1 = rc*dx[0] - 0.5*dx[0];
		volume = 4.0*PI*(r2 - r1)*(r2*r2 + r1*r2 + r1*r1)/3.0;
	}
	return volume;
}



