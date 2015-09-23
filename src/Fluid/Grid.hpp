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

#include <array>
#include <utility>
#include <memory>

#include "Torch/Common.hpp"
#include "Torch/Constants.hpp"
#include "Torch/Parameters.hpp"
#include "GridCell.hpp"
#include "PartitionManager.hpp"
#include "GridCellCollection.hpp"

class Boundary;

class Bound {
public:
	int face;
	Condition condition;
	int targetProcessor;
	std::vector<int> ghostCellIDs;


	Bound(int face, const Condition bcond, int target_proc = 0);
};

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
public:
	std::vector<int> m_causalIndices;
	std::vector<Bound> m_boundaries;
	std::array<int, 3> ncells = std::array<int, 3>{ 1, 1, 1 };
	std::array<int, 3> coreCells = std::array<int, 3>{ 1, 1, 1 };
	std::array<double, 3> dx = std::array<double, 3>{ 0, 0, 0 };
	double sideLength = 0;
	int spatialOrder = 0;
	Geometry geometry = Geometry::CARTESIAN;
	double currentTime = 0; //!< current time in simulation.
	double deltatime = 0; //!< the time to increment. TO-DO: move to Integrator.


	Grid();

	void initialise(std::shared_ptr<Constants> consts, const GridParameters& gp);
	void applyBCs();

	GridCell& getCell(int id);
	const GridCell& getCell(int id) const;
	Looper getIterable(const std::string& name);
	ConstLooper getIterable(const std::string& name) const;
	const std::vector<GridCell>& getCells() const;
	const std::vector<GridJoin>& getJoins(int dim) const;
	std::vector<GridCell>& getCells();
	std::vector<GridJoin>& getJoins(int dim);
	PartitionManager& getPartitionManager();
	std::vector<int>& getCausalIndices();
	std::vector<int>& getOrderedIndices(const std::string& name);

	std::vector<Bound>& getBoundaries();

	void addOrderedIndex(const std::string& name, int index);

	int getLeftX() const;
	int getRightX() const;
	void setLeftX(int x);
	void setRightX(int x);

	bool cellExists(int id) const;

	// Traversal.
	int left(int dim, int fromCellID);
	int right(int dim, int fromCellID);
	int ghostLeft(int dim, int fromCellID);
	int ghostRight(int dim, int fromCellID);
	GridCell& left(int dim, GridCell& fromCell);
	GridCell& right(int dim, GridCell& fromCell);
	GridCell& ghostLeft(int dim, GridCell& fromCell);
	GridCell& ghostRight(int dim, GridCell& fromCell);
	GridCell& joinedLeft(GridJoin& join);
	GridCell& joinedRight(GridJoin& join);
	GridJoin& leftJoin(int dim, GridCell& fromCell);
	GridJoin& rightJoin(int dim, GridCell& fromCell);

	// Traversal.
	int locate(int cx, int cy, int cz);
	Coords nearestCoord(const Coords& original);
	bool joinExists(int dim, std::array<int, 3>& joinIDs);
	bool withinGrid(const Coords& coords);
	Coords nextSnakeCoords(const Coords& fromCoords, const Coords& sourceCoords, const int dxc, const int dyc, const int dyz);
	Coords nextCausalCoords(const Coords& fromCoords, const Coords& sourceCoords);
	int traverse1D(const int dim, const int dc, int fromCellID);
	int traverse3D(const int d1, const int d2, const int d3, const int dc1, const int dc2, const int dc3, int fromCellID);
	int traverseOverJoins1D(const int dim, const int dc, int fromCellID);
	int traverseOverJoins3D(const int d1, const int d2, const int d3, const int dc1, const int dc2, const int dc3, int fromCellID);
	int nextCell2D(const int plane, int fromCellID);
	int nextSnake(int fromCellID, int sourceCellID, const int dxc, const int dyc, const int dyz, int nd);
	int nextCausal(int fromCellID, int sourceCellID, int nd);

	// Coord transform.
	int flatIndex(int ci, int cj, int ck);
	Coords unflatCoords(int flat_index);

	// Creation.
	void weakLink(const int dim, int lcellID, int rcellID);
	void link(const int dim, int lcellID, int rcellID);
	void boundaryLink(Bound& boundary);
	void boundaryLinkDeeper(Bound& boundary);
	void buildCells();
	void buildCausal(const Coords& sourceCoords);
	void buildBoundaries(const std::array<Condition, 3>& leftBC, const std::array<Condition, 3>& rightBC);
	double computeCellVolume(double rc, const Vec3& dx, Geometry geometry, int nd);
	double computeJoinArea(const Vec3& xj, const int dim, const Vec3& dx, Geometry geometry, int nd);
	int getRayPlane(const Vec3& xc, const Vec3& xs) const;
	void calculateNearestNeighbours(const std::array<double, 3>& star_pos);
private:
	std::shared_ptr<Constants> m_consts = nullptr;
	PartitionManager partition;
	GridCellCollection m_cellCollection;
	std::vector<GridCell>& m_cells = m_cellCollection.getCellVector();
	std::map<std::string, std::vector<int>> orderedIndices;
	std::array<std::vector<GridJoin>, 3> m_joins = std::array<std::vector<GridJoin>, 3>{ std::vector<GridJoin>(), std::vector<GridJoin>(), std::vector<GridJoin>() };
	int m_leftX = 0;
	int m_rightX = 0;
};


#endif // GRID_HPP_
