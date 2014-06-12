/**
 * Provides the Grid3D class.
 * @file grid3d.hpp
 *
 * @author Harrison Steggles
 *
 * @date 13/01/2014 - the first version.
 * @date 16/01/2014 - removed old Boundary class. New one now holds all ghostcells associated with a grid face. Two types of Boundary can be
 * attached to a simulation grid: ExternalBoundary and Partition.
 * @date 29/01/2014 - multiple processor capabilities added.
 * @date 04/02/2014 - added pointer to first GridCell in causal iteration and a method (buildCausal) for setting up this loop.
 * @date 04/02/2014 - added pointers to last GridJoins and GridCells in "next" lists so that adding GridCells to end of list is a LOT faster.
 * This has removed a major pre-simulation bottleneck.
 * @date 04/01/2014 - arguments now passed by const reference when appropriate.
 * @date 05/02/2014 - fixed nextSnake bug causing infinite loop in 3D buildCausal.
 * @date 11/02/2014 - fixed bug in addGridCell causing 1D (instead of 3D) grids to be set up when using mpi with more than 1 core.
 * @date 21/02/2014 - added InputOutput::progressBar to loop that builds grid and fixed bug in traverseOverJoins3D incorrectly returning
 * null.
 * @date 24/02/2014 - Grid3D now decides to override user specified boundary conditions when particular geometries are set up.
 * @date 14/04/2014 - Grid3D parameters now all held within a GridParameters object.
 */

#ifndef GRID3D_H
#define GRID3D_H

#include "constants.hpp"

#include <vector>

class GridCell;
class GridJoin;
class Boundary;
class GridParameters;
class Scalings;
class MPIHandler;
/** @class Grid3D
 * @brief The Grid3D class holds the state information of the fluid and provides methods for traversing and locating {@link GridCell}s.
 *
 * A Grid3D object has access to a doubly linkedlist of GridCell objects that hold fluid and radiation
 * state information. It also links the grid of simulation GridCell objects with appropriate Boundary
 * derived objects. Traversal can either be plough-like or causal wrt a Star object.
 *
 * The grid coordinates are integers that start at 0. The coordinate in simulation units is given by (xc[dim] + 0.5)*dx[dim].
 *
 * @author      Harrison Steggles
 * @version     0.5, 24/01/2014
 */
class Grid3D
{
public:
	static int s_total; //!< Number of Grid3D objects in the program.
	GridParameters* gparams; //!< Pointer to object containing grid parameters.
	Scalings* scale; //!< Pointer to object containing scalings between computational and physical units.
	GridCell* fcell; //!< Pointer to first GridCell object in simulation grid.
	GridCell* lcell; //!< Pointer to last GridCell object in simulation grid.
	GridJoin* fjoin[3]; //!< Pointers to first GridJoin object in linked list in each dimension.
	GridJoin* ljoin[3]; //!< Pointers to last GridJoin object in linked list in each dimension.
	std::vector<Boundary*> boundaries; //!< Pointers to the Boundary objects connected to Grid3D.
	double dx[3]; //!< Cell width in each dimension.
	double currentTime; //!< current time in simulation.
	double deltatime; //!< the time to increment. TO-DO: move to Integrator.

	//Structors.
	Grid3D(const int spatialOrder, const GridParameters& gp, const Scalings& sc, MPIHandler& mpih);
	~Grid3D();

	//Structure linking methods.
	void buildCausal(GridCell* fcausal);

	//Traversal methods.
	GridCell* traverse1D(const int& dim, const int& dc, GridCell* cptr);
	GridCell* traverse3D(const int& d1, const int& d2, const int& d3, const int& dc1, const int& dc2, const int& dc3, GridCell* cptr);
	GridCell* traverseOverJoins1D(const int& dim, const int& dc, GridCell* cptr);
	GridCell* traverseOverJoins3D(const int& d1, const int& d2, const int& d3, const int& dc1, const int& dc2, const int& dc3, GridCell* cptr);
	GridCell* locate(const int& x, const int& y, const int& z);
	GridCell* nextCell2D(const int& plane, GridCell* cptr);
	GridCell* nextCell3D(GridCell* cptr);
	GridCell* nextSnake(GridCell* cptr, GridCell* srcptr, const int& dxc, const int& dyc, const int& dyz);
	GridCell* nextCausal(GridCell* cptr, GridCell* srcptr);

private:
	//Structure linking methods.
	void link(const int& dim, GridCell* lcptr, GridCell* rcptr);
	void weakLink(const int& dim, GridCell* lcptr, GridCell* rcptr);
	void boundaryLink(Boundary* bptr);

	//Add methods.
	void addGridCell(GridCell* newcptr);
	void addJoinToList(GridJoin* jptr, const int& dim);
	void addGridCellToList(GridCell* cptr);

	//Geometry calculation methods.
	void vol_cell(GridCell* cptr);
	double area_join(double xj[], const int& dim);
};

#endif
