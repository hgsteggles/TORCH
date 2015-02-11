/**
 * Provides the Grid class.
 * @file Grid.h
 *
 * @author Harrison Steggles
 *
 * @date 13/01/2014 - the first version.
 * @date 16/01/2014 - removed old Boundary class. New one now holds all ghostcells associated with a grid face. Two types of Boundary can be
 * attached to a simulation grid: ExternalBoundary and Partition.
 * @date 29/01/2014 - multiple processor capabilities added.
 * @date 04/02/2014 - added pointer to first GridCell in causal iteration and a method (buildCausal) for setting up this loop.
 * @date 04/02/2014 - added pointers to last GridJoins and GridCells in "next" lists so that adding GridCells to end of list is a lot faster.
 * This has removed a major pre-simulation bottleneck.
 * @date 04/01/2014 - arguments now passed by const reference when appropriate.
 * @date 05/02/2014 - fixed nextSnake bug causing infinite loop in 3D buildCausal.
 * @date 11/02/2014 - fixed bug in addGridCell causing 1D (instead of 3D) grids to be set up when using mpi with more than 1 core.
 * @date 21/02/2014 - added InputOutput::progressBar to loop that builds grid and fixed bug in traverseOverJoins3D incorrectly returning
 * null.
 * @date 24/02/2014 - Grid3D now decides to override user specified boundary conditions when particular geometries are set up.
 * @date 14/04/2014 - Grid3D parameters now all held within a GridParameters object.
 * @date 21/07/2014 - replaced const reference to int/double with a const copy.
 *
 * @date 24/11/2014 - new class Grid is Grid3D without the creation. The creation is now handled separately by GridFactory. A Grid now
 * has a "DualIntrusiveContainer" of GridCells so that two separate lists can be traversed and the entire list can be traversed (useful
 * for separating stellar wind cells from the rest).
 */

#ifndef GRID_HPP_
#define GRID_HPP_

#include "Common.hpp"
#include "Constants.hpp"
#include "GridCell.hpp"

class Boundary;

/**
 * @class Grid
 *
 * @brief The Grid class holds all the GridCells and GridJoins along with Boundary's that enclose the simulation domain
 * and provides methods for traversing and locating {@link GridCell}s.
 *
 * A Grid contains a doubly linked list of GridCells that hold fluid and radiation state information.
 * Traversal can either be plough-like, causal with respect to a Star, or along each dimension starting from a specified GridCell.
 *
 * The grid coordinates are in units of the GridCell widths and are located at the centre of each GridCell. At the origin (corner
 * of the GridCell) the coordinates are {0, 0, 0}. The coordinate in simulation units is given by xc[dim]*dx[dim]. Use a Converter
 * to convert to cgs units.
 *
 * @version 0.8, 24/11/2014
 */
class Grid {
	using GridCellPair = std::pair<GridCell* const, GridCell* const>;
public:
	DualIntrusiveContainer<GridCell> m_cells;
	std::array<JoinContainer, 3> m_joins;
	std::array<Boundary*, 3> m_leftBoundaries = std::array<Boundary*, 3>{ nullptr, nullptr, nullptr };
	std::array<Boundary*, 3> m_rightBoundaries = std::array<Boundary*, 3>{ nullptr, nullptr, nullptr };

	std::array<int, 3> ncells = std::array<int, 3>{ 1, 1, 1 };
	std::array<int, 3> coreCells = std::array<int, 3>{ 1, 1, 1 };
	double sideLength = 0;
	int spatialOrder = 0;
	Geometry geometry = Geometry::CARTESIAN;
	std::array<double, 3> dx = std::array<double, 3>{ 0, 0, 0 };

	double currentTime = 0; //!< current time in simulation.
	double deltatime = 0; //!< the time to increment. TO-DO: move to Integrator.


	Grid();
    Grid(Grid&&);
    Grid& operator=(Grid&&);
	~Grid();

	int getLeftX() const;
	int getRightX() const;
	void setLeftX(int x);
	void setRightX(int x);

	const CellContainer& getCells() const;
	const CellContainer& getWindCells() const;
	const CellContainer& getCausalCells() const;
	const JoinContainer& getJoins(int dim) const;
	const std::array<Boundary*, 3>& getLeftBoundaries() const;
	const std::array<Boundary*, 3>& getRightBoundaries() const;

	CellContainer& getCells();
	CellContainer& getWindCells();
	CellContainer& getCausalCells();
	JoinContainer& getJoins(int dim);
	std::array<Boundary*, 3>& getLeftBoundaries();
	std::array<Boundary*, 3>& getRightBoundaries();

	// Traversal.
	static GridCell* traverse1D(const int dim, const int dc, GridCell* const cptr);
	static GridCell* traverse3D(const int d1, const int d2, const int d3, const int dc1, const int dc2, const int dc3, GridCell* const cptr);
	static GridCell* traverseOverJoins1D(const int dim, const int dc, GridCell* const cptr);
	static GridCell* traverseOverJoins3D(const int d1, const int d2, const int d3, const int dc1, const int dc2, const int dc3, GridCell* const cptr);
	static GridCell* nextCell2D(const int plane, GridCell* const cptr);
	static GridCell* nextSnake(GridCell* const cptr, GridCell* srcptr, const int dxc, const int dyc, const int dyz, int nd);
	static GridCell* nextCausal(GridCell* const cptr, GridCell* srcptr, int nd);
	static GridCell* locate(const Coords& x, GridCell* const fromCell);
	static GridCell* getNearestCell(const Coords& x, const GridCellPair& firstLast);
private:
	int m_leftX = 0;
	int m_rightX = 0;

    Grid(const Grid&);
    Grid& operator=(const Grid&);
};



#endif // GRID_HPP_
