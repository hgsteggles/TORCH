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
 */

#ifndef GRID3D_H
#define GRID3D_H

#include "constants.hpp"

#include <vector>

class GridCell;
class GridJoin;
class Boundary;
class GridParameters;
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
	static int s_total; //!< Number of Grid3D objects in program.
	GridCell* fcell; //!< Pointer to first GridCell object in simulation grid.
	GridCell* lcell; //!< Pointer to last GridCell object in simulation grid.
	GridJoin* fjoin[3]; //!< Pointers to first GridJoin object in linked list in each dimension.
	GridJoin* ljoin[3];
	double dx[3]; //!< Cell width in each dimension.
	Geometry GEOMETRY; //!< Enum which notes what coordinate system is in use.
	int ND, NU, NR, TOTNCELLS[3], NCELLS[3], ORDER_S, ORDER_T;
	std::vector<Boundary*> boundaries;
	double currentTime;
	double deltatime;
	/**
	 * @brief Grid3D constructor.
	 * Sets Grid3D parameters equal to parameters owned by GridParameters. Uses parameters to construct
	 * a 3D double-linked list of GridCell objects each of which is linked to GridJoin objects positioned
	 * on the GridCell faces. ExternalBoundary and Partition objects are instantiated and held in a vector
	 * of Boundary pointers. The Boundary objects are linked to the faces of the grid via GridJoin objects
	 * between GridCell objects.
	 *
	 * @param gp
	 * @param mpih
	 */
	Grid3D(const GridParameters& gp, MPIHandler& mpih);
	~Grid3D();
	/**
	 * @brief Double-links the passed GridCell objects to each other and a GridJoin object.
	 * Links together two GridCell objects via GridCell::left and GridCell::right and also links them
	 * both to a common GridJoin object via GridCell::rjoin and GridCell::ljoin. GridJoin contains links
	 * to both GridCell objects via GridJoin::rcell and GridJoin::lcell. The coordinates of the
	 * GridJoin object is set according to those of the GridCell objects.
	 *
	 * @param dim Dimension across which to set up links.
	 * @param lcptr Pointer to left GridCell object.
	 * @param rcptr Pointer to right GridCell object.
	 */
	void link(const int& dim, GridCell* lcptr, GridCell* rcptr);
	/**
	 * @brief Double-links the passed GridCell objects to a GridJoin object.
	 * Links two GridCell objects to a common GridJoin object via GridCell::rjoin and GridCell::ljoin.
	 * GridJoin contains links to both GridCell objects via GridJoin::rcell and GridJoin::lcell. The
	 * coordinates of the GridJoin object is set according to those of the GridCell objects.
	 * @param dim Dimension across which to set up links.
	 * @param lcptr Pointer to left GridCell object.
	 * @param rcptr Pointer to right GridCell object.
	 */
	void weakLink(const int& dim, GridCell* lcptr, GridCell* rcptr);
	/**
	 * @brief Links GridCell objects on a grid face with GridCell objects in a Boundary object via a
	 * GridJoin object. Each GridCell object on a grid face is linked via a GridJoin to a GridCell object
	 * owned by a Boundary. The GridCell objects are not double-linked to each other because traversal
	 * over the simulation grid is achieved by detecting NULL pointers on the boundary.
	 * @param bptr Pointer to the Boundary object to be linked.
	 */
	void boundaryLink(Boundary* bptr);
	/**
	 * @brief Adds a GridCell to the 3D double-linked list of GridCell objects.
	 * @param newcptr Pointer to the GridCell object to be added.
	 */
	void addGridCell(GridCell* newcptr);
	/**
	 * @brief Adds a GridJoin object to a linked list owned by Grid3D.
	 * Useful for iterating over the objects and cleaning them up in the Grid3D destructor.
	 * @param jptr Pointer to GridJoin object to add.
	 * @param dim Coordinate axis that is normal to the GridJoin object.
	 */
	void addJoinToList(GridJoin* jptr, const int& dim);
	/**
	 * @brief Adds a GridCell object to a linked list owned by Grid3D.
	 * Useful for iterating over the objects and cleaning them up in the Grid3D destructor.
	 * @param cptr Pointer to GridCell object to add.
	 */
	void addGridCellToList(GridCell* cptr);
	/**
	 * @brief Traverses over GridCell objects via GridCell::right and GridCell::left in 1D.
	 * @param dim The direction to traverse (0, 1, 2).
	 * @param dc No. of GridCell objects to traverse over. Passing 1 would return the adjacent GridCell object.
	 * @param cptr Pointer to the GridCell object to start from.
	 * @return The GridCell object traversed to. Will return <code> NULL </code> if it doesn't exist.
	 */
	GridCell* traverse1D(const int& dim, const int& dc, GridCell* cptr);
	/**
	 * @brief Traverses over GridCell objects via GridCell::right and GridCell::left in 3D.
	 * @param d1 Direction 1.
	 * @param d2 Direction 2.
	 * @param d3 Direction 3.
	 * @param dc1 No. of GridCell objects to traverse along direction d1.
	 * @param dc2 No. of GridCell objects to traverse along direction d2.
	 * @param dc3 No. of GridCell objects to traverse along direction d3.
	 * @param cptr Pointer to the GridCell object to start from.
	 * @return The GridCell object traverse to. Will return <code> NULL </code> if it doesn't exist.
	 */
	GridCell* traverse3D(const int& d1, const int& d2, const int& d3, const int& dc1, const int& dc2, const int& dc3, GridCell* cptr);
	/**
	 * @brief Traverses over GridCell objects via GridCell::rjoin and GridCell::ljoin in 1D.
	 * @param dim The direction to traverse (0, 1, 2).
	 * @param dc No. of GridCell objects to traverse over. Passing 1 would return the adjacent GridCell object.
	 * @param cptr Pointer to the GridCell object to start from.
	 * @return The GridCell object traversed to. Will return <code> NULL </code> if it doesn't exist.
	 */
	GridCell* traverseOverJoins1D(const int& dim, const int& dc, GridCell* cptr);
	/**
	 * @brief Traverses over GridCell objects via GridCell::rjoin and GridCell::ljoin in 1D.
	 * @param d1 Direction 1.
	 * @param d2 Direction 2.
	 * @param d3 Direction 3.
	 * @param dc1 No. of GridCell objects to traverse along direction d1.
	 * @param dc2 No. of GridCell objects to traverse along direction d2.
	 * @param dc3 No. of GridCell objects to traverse along direction d3.
	 * @param cptr Pointer to the GridCell object to start from.
	 * @return The GridCell object traverse to. Will return <code> NULL </code> if it doesn't exist.
	 */
	GridCell* traverseOverJoins3D(const int& d1, const int& d2, const int& d3, const int& dc1, const int& dc2, const int& dc3, GridCell* cptr);
	/**
	 * @brief Locates GridCell object at the specified grid coordinates.
	 * @param x x grid coordinate.
	 * @param y y grid coordinate.
	 * @param z z grid coordinate.
	 * @return The located GridCell object if it exists. If not then <code> NULL </code>.
	 */
	GridCell* locate(const int& x, const int& y, const int& z);
	/**
	 * @brief Provides the next GridCell object along a path that traverses a 2D plane which includes the passed GridCell object.
	 * @param plane The 2D plane normal.
	 * @param cptr Pointer to previous GridCell object
	 * @return Pointer to next GridCell object. If there are no more GridCell objects on this path then <code> NULL </code> is returned.
	 */
	GridCell* nextCell2D(const int& plane, GridCell* cptr);
	/**
	 * @brief Provides the next GridCell object along a path that traverses the entire grid.
	 * @param cptr Pointer to previous GridCell object
	 * @return Pointer to next GridCell object. If there are no more GridCell objects on this path then <code> NULL </code> is returned.
	 */
	GridCell* nextCell3D(GridCell* cptr);
	/**
	 * @brief Provides next GridCell object along a specific path.
	 * The path starts at the GridCell object pointed to by srcptr and travels along the positive or negative x direction (sign of dxc)
	 * until a NULL pointer is encountered. The next GridCell object in this case would be along the y direction (specified by dyc) from
	 * a GridCell that has the x and z coordinates of the GridCell pointed to by srcptr and the y coordinate of the previous GridCell
	 * object. At some point traversing in the y direction will encounter a NULL pointer a so step is taken in the z direction starting
	 * from a GridCell object with x and y coordinates of the GridCell pointed to by srcptr and z coordinate of the previous GridCell
	 * object. E.g. (dxc,dyc,dzc)=(1,1,1) leads to all GridCell objects with x,y,z coordinate equal to or greater than the GridCell
	 * pointed to by srcptr being sampled along this route.
	 *
	 * @param srcptr The GridCell object that traversal in all directions starts from.
	 * @param cptr Pointer to previous GridCell object.
	 * @param dxc GridCell objects to traverse in the x direction at every step.
	 * @param dyc GridCell objects to traverse in the y direction at every step.
	 * @param dyz GridCell objects to traverse in the z direction at every step.
	 * @return Pointer to next GridCell object. If there are no more GridCell objects on this path then <code> NULL </code> is returned.
	 */
	GridCell* nextSnake(GridCell* cptr, GridCell* srcptr, const int& dxc, const int& dyc, const int& dyz);
	/**
	 * @brief Provides next GridCell object along a causal path from a specified GridCell object.
	 * @param cptr Pointer to previous GridCell object.
	 * @param srcptr Pointer to the first GridCell along the path.
	 * @return Pointer to next GridCell object. If there are no more GridCell objects on this path then <code> NULL </code> is returned.
	 */
	GridCell* nextCausal(GridCell* cptr, GridCell* srcptr);
	void buildCausal(GridCell* fcausal);
	/**
	 * @brief Calculates the area between 2 GridCell objects given its location and the geometry of Grid3D.
	 * @param xj Coordinates of the GridJoin between the 2 GridCell objects.
	 * @param dim Direction normal to the face of the GridJoin.
	 * @return The area.
	 */
	double area_join(double xj[], const int& dim);
	/**
	 * @brief Calculates GridCell::vol given its location and  the geometry of Grid3D.
	 * @param cptr The GridCell object for which the volume will be calculated.
	 */
	void vol_cell(GridCell* cptr);
};

#endif
