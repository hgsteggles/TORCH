/**
 * @file grid3d.cpp
 */

#include "grid3d.h"
#include "parameters.h"
#include "mpihandler.h"
#include "gridcell.h"
#include "external.h"
#include "partition.h"
#include "io.h"

#include <cstring>
#include <cmath>

/**
 * @brief Grid3D constructor.
 * Sets Grid3D parameters equal to parameters owned by GridParameters. Uses parameters to construct
 * a 3D double-linked list of GridCell objects each of which is linked to GridJoin objects positioned
 * on the GridCell faces. ExternalBoundary and Partition objects are instantiated and held in a vector
 * of Boundary pointers. The Boundary objects are linked to the faces of the grid via GridJoin objects
 * between GridCell objects.
 *
 * @param spatialOrder
 * @param gp
 * @param sc
 * @param mpih
 */
Grid3D::Grid3D(const int spatialOrder, const GridParameters& gp, MPIHandler& mpih) : fcell(NULL), lcell(NULL), currentTime(0.0), deltatime(0.0) {
	gparams = new GridParameters(gp);
	gparams->NCELLS[1] = gparams->ND > 1 ? gparams->NCELLS[1] : 1;
	gparams->NCELLS[2] = gparams->ND > 2 ? gparams->NCELLS[2] : 1;
	std::memcpy( (void*)gparams->CORECELLS, (void*)gparams->NCELLS, 3 * sizeof(int) );
	if (gparams->NCELLS[0]%mpih.nProcessors() == 0)
		gparams->CORECELLS[0] = gparams->NCELLS[0] / mpih.nProcessors();
	else {
		std::cout << "ERROR: No. of GridCells (" << gparams->NCELLS << ") in x dir not divisible by no. of processors (" << mpih.nProcessors() << ")." << '\n';
		exit(EXIT_FAILURE);
	}
	for(int i = 0; i < 3; ++i){
		fjoin[i] = NULL;
		ljoin[i] = NULL;
	}
	dx[0] = gparams->SIDE_LENGTH/(double)gparams->NCELLS[0];
	dx[1] = gparams->SIDE_LENGTH/(double)gparams->NCELLS[1];
	dx[2] = gparams->SIDE_LENGTH/(double)gparams->NCELLS[2];

	InputOutput::initProgressBar("Building Grid", mpih);
	for(int i = 0; i < (gparams->CORECELLS[0]*gparams->CORECELLS[1]*gparams->CORECELLS[2]); ++i){
		InputOutput::progressBar((100*i/(double)(gparams->CORECELLS[0]*gparams->CORECELLS[1]*gparams->CORECELLS[2])), 5, mpih);
		GridCell* newcptr = new GridCell();
		if (i == 0)
			newcptr->set_xcs(mpih.getRank()*gparams->CORECELLS[0], 0, 0);
		addGridCell(newcptr);
		vol_cell(newcptr);
	}
	InputOutput::endProgressBar(mpih);
	std::cout << mpih.cname() <<  "GRID3D: set up with " << gparams->CORECELLS[0]*gparams->CORECELLS[1]*gparams->CORECELLS[2] << " cells." << '\n';
	if(gparams->ND > 0) {
		if (mpih.getRank() == 0) {
			if (gparams->GEOMETRY == SPHERICAL || gparams->GEOMETRY == CYLINDRICAL)
				boundaries.push_back(new ExternalBoundary(0, spatialOrder+1, REFLECTING, this));
			else
				boundaries.push_back(new ExternalBoundary(0, spatialOrder+1, gp.LBCondition[0], this));
		}
		else {
			int destination = mpih.getRank()-1;
			boundaries.push_back(new Partition(0, spatialOrder+1, this, destination, mpih));
		}
		boundaryLink(boundaries[0]);
		if (mpih.getRank() == mpih.nProcessors()-1) {
			boundaries.push_back(new ExternalBoundary(3, spatialOrder+1, gp.RBCondition[0], this));
		}
		else {
			boundaries.push_back(new Partition(3, spatialOrder+1, this, mpih.getRank()+1, mpih));
		}
		boundaryLink(boundaries[1]);
	}
	if(gparams->ND > 1) {
		boundaries.push_back(new ExternalBoundary(1, spatialOrder+1, gp.LBCondition[1], this));
		boundaryLink(boundaries[2]);
		boundaries.push_back(new ExternalBoundary(4, spatialOrder+1, gp.RBCondition[1], this));
		boundaryLink(boundaries[3]);
	}
	if(gparams->ND > 2) {
		boundaries.push_back(new ExternalBoundary(2, spatialOrder+1, gp.LBCondition[2], this));
		boundaryLink(boundaries[4]);
		boundaries.push_back(new ExternalBoundary(5, spatialOrder+1, gp.RBCondition[2], this));
		boundaryLink(boundaries[5]);
	}
	//boundaries[1]->ghostcells[0][62]->printInfo();

	s_total++;
}
Grid3D::~Grid3D(){ 
	GridJoin* nextjptr = NULL;
	GridCell* nextcptr = NULL;
	for (int i = 0; i < (int)boundaries.size(); ++i) {
		delete boundaries[i];
		boundaries[i] = NULL;
	}
	for(int dim = 0; dim < gparams->ND; ++dim){
		for(GridJoin* jptr = fjoin[dim]; jptr != NULL; jptr = nextjptr){
			nextjptr = jptr->next;
			delete jptr;
		}
	}
	for(GridCell* cptr = fcell; cptr != NULL; cptr = nextcptr){
		nextcptr = cptr->next;
		delete cptr;
	}
}

int Grid3D::s_total = 0;

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
void Grid3D::link(const int& dim, GridCell* lcptr, GridCell* rcptr){
	if(lcptr != NULL && rcptr != NULL){
		rcptr->left[dim] = lcptr;
		lcptr->right[dim] = rcptr;
		GridJoin* jptr = new GridJoin();
		rcptr->ljoin[dim] = jptr;
		lcptr->rjoin[dim] = jptr;
		jptr->rcell = rcptr;
		jptr->lcell = lcptr;
		for(int i = 0; i < 3; ++i)
			jptr->xj[i] = jptr->rcell->xc[i];
		jptr->xj[dim] = jptr->rcell->xc[dim] - 0.5;
		jptr->area = area_join(jptr->xj, dim);
		addJoinToList(jptr, dim);
	}
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
void Grid3D::weakLink(const int& dim, GridCell* lcptr, GridCell* rcptr){
	if(lcptr != NULL && rcptr != NULL){
		GridJoin* jptr = new GridJoin();
		rcptr->ljoin[dim] = jptr;
		lcptr->rjoin[dim] = jptr;
		jptr->rcell = rcptr;
		jptr->lcell = lcptr;
		for(int i = 0; i < 3; ++i)
			jptr->xj[i] = jptr->rcell->xc[i];
		jptr->xj[dim] = jptr->rcell->xc[dim] - 0.5;
		jptr->area = area_join(jptr->xj, dim);
		addJoinToList(jptr, dim);
	}
}

/**
 * @brief Links GridCell objects on a grid face with GridCell objects in a Boundary object via a
 * GridJoin object. Each GridCell object on a grid face is linked via a GridJoin to a GridCell object
 * owned by a Boundary. The GridCell objects are not double-linked to each other because traversal
 * over the simulation grid is achieved by detecting NULL pointers on the boundary.
 * @param bptr
 * Pointer to the Boundary object to be linked.
 */
void Grid3D::boundaryLink(Boundary* bptr) {
	int dim = bptr->face%3;
	if (bptr->face < 3) {
		for (GridCell* cptr = fcell; cptr != NULL; cptr = nextCell2D(dim, cptr)) {
			GridCell* bcptr = NULL;
			if (dim == 0)
				bcptr = bptr->ghostcells[cptr->xc[1]][cptr->xc[2]];
			else if (dim == 1)
				bcptr = bptr->ghostcells[cptr->xc[0]-fcell->xc[0]][cptr->xc[2]];
			else
				bcptr = bptr->ghostcells[cptr->xc[0]-fcell->xc[0]][cptr->xc[1]];
			for (int i = 0; i < 3; ++i)
				bcptr->xc[i] = cptr->xc[i];
			bcptr->xc[dim] = cptr->xc[dim] - 1;
			weakLink(dim, bcptr, cptr);
		}
	}
	else {
		for (GridCell* cptr = traverse1D(dim, gparams->CORECELLS[dim]-1, fcell); cptr != NULL; cptr = nextCell2D(dim, cptr)) {
			GridCell* bcptr = NULL;
			if (dim == 0)
				bcptr = bptr->ghostcells[cptr->xc[1]][cptr->xc[2]];
			else if (dim == 1)
				bcptr = bptr->ghostcells[cptr->xc[0]-fcell->xc[0]][cptr->xc[2]];
			else
				bcptr = bptr->ghostcells[cptr->xc[0]-fcell->xc[0]][cptr->xc[1]];
			for (int i = 0; i < 3; ++i)
				bcptr->xc[i] = cptr->xc[i];
			bcptr->xc[dim] = cptr->xc[dim] + 1;
			weakLink(dim, cptr, bcptr);
		}
	}
}

/**
 * @brief Traverses over GridCell objects via GridCell::right and GridCell::left in 1D.
 * @param dim
 * The direction to traverse (0, 1, 2).
 * @param dc
 * No. of GridCell objects to traverse over. Passing 1 would return the adjacent GridCell object.
 * @param cptr
 * Pointer to the GridCell object to start from.
 * @return The GridCell object traversed to. Will return <code> NULL </code> if it doesn't exist.
 */
GridCell* Grid3D::traverse1D(const int& dim, const int& dxc, GridCell* cptr){
	GridCell* newcptr = cptr;
	for(int i = 0; i < abs(dxc) && newcptr != NULL; ++i){
		if(dxc > 0)
			newcptr = newcptr->right[dim];
		else
			newcptr = newcptr->left[dim]; 
	}
	return newcptr;
}

/**
 * @brief Traverses over GridCell objects via GridCell::right and GridCell::left in 3D.
 * blach
 * @param d1
 * Direction 1.
 * @param d2
 * Direction 2.
 * @param d3
 * Direction 3.
 * @param dc1
 * No. of GridCell objects to traverse along direction d1.
 * @param dc2
 * No. of GridCell objects to traverse along direction d2.
 * @param dc3
 * No. of GridCell objects to traverse along direction d3.
 * @param cptr
 * Pointer to the GridCell object to start from.
 * @return The GridCell object traverse to. Will return <code> NULL </code> if it doesn't exist.
 */
GridCell* Grid3D::traverse3D(const int& d1, const int& d2, const int& d3, const int& dc1, const int& dc2, const int& dc3, GridCell* cptr){
	return traverse1D(d3, dc3, traverse1D(d2, dc2, traverse1D(d1, dc1, cptr)));
}

/**
 * @brief Traverses over GridCell objects via GridCell::rjoin and GridCell::ljoin in 1D.
 * @param dim
 * The direction to traverse (0, 1, 2).
 * @param dc
 * No. of GridCell objects to traverse over. Passing 1 would return the adjacent GridCell object.
 * @param cptr
 * Pointer to the GridCell object to start from.
 * @return The GridCell object traversed to. Will return <code> NULL </code> if it doesn't exist.
 */
GridCell* Grid3D::traverseOverJoins1D(const int& dim, const int& dxc, GridCell* cptr){
	GridCell* newcptr = cptr;
	for(int i = 0; i < abs(dxc) && newcptr != NULL; ++i){
		if (dxc > 0 && newcptr->rjoin[dim] != NULL)
			newcptr = newcptr->rjoin[dim]->rcell;
		else if (dxc < 0 && newcptr->ljoin[dim] != NULL)
			newcptr = newcptr->ljoin[dim]->lcell;
		else
			newcptr = NULL;
	}
	return newcptr;
}

/**
 * @brief Traverses over GridCell objects via GridCell::rjoin and GridCell::ljoin in 1D.
 * @param d1
 * Direction 1.
 * @param d2
 * Direction 2.
 * @param d3
 * Direction 3.
 * @param dc1
 * No. of GridCell objects to traverse along direction d1.
 * @param dc2
 * No. of GridCell objects to traverse along direction d2.
 * @param dc3
 * No. of GridCell objects to traverse along direction d3.
 * @param cptr
 * Pointer to the GridCell object to start from.
 * @return The GridCell object traverse to. Will return <code> NULL </code> if it doesn't exist.
 */
GridCell* Grid3D::traverseOverJoins3D(const int& d1, const int& d2, const int& d3, const int& dc1, const int& dc2, const int& dc3, GridCell* cptr){
	GridCell* newcptr = traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d1, dc1, cptr)));
	if (newcptr == NULL)
		newcptr = traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d3, dc3, cptr)));
	if (newcptr == NULL)
		newcptr = traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d2, dc2, cptr)));
	if (newcptr == NULL)
		newcptr = traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d3, dc3, cptr)));
	if (newcptr == NULL)
		newcptr = traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d1, dc1, cptr)));
	if (newcptr == NULL)
		newcptr = traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d2, dc2, cptr)));
	return newcptr;
}

/**
 * @brief Locates GridCell object at the specified grid coordinates.
 * @param x
 * x grid coordinate.
 * @param y
 * y grid coordinate.
 * @param z
 * z grid coordinate.
 * @return The located GridCell object if it exists. If not then <code> NULL </code>.
 */
GridCell* Grid3D::locate(const int& x, const int& y, const int& z){
	GridCell* cptr = fcell;
	for(int i = 0; (cptr != NULL) && (cptr->xc[0] != x); ++i)
		cptr = cptr->right[0];
	for(int i = 0; (cptr != NULL) && (cptr->xc[1] != y); ++i)
		cptr = cptr->right[1];
	for(int i = 0; (cptr != NULL) && (cptr->xc[2] != z); ++i)
		cptr = cptr->right[2];
	return cptr;
}

/**
 * @brief Adds a GridJoin object to a linked list owned by Grid3D.
 * Useful for iterating over the objects and cleaning them up in the Grid3D destructor.
 * @param jptr
 * Pointer to GridJoin object to add.
 * @param dim
 * Coordinate axis that is normal to the GridJoin object.
 */
void Grid3D::addJoinToList(GridJoin* jptr, const int& dim){
	if(fjoin[dim] != NULL){
		ljoin[dim]->next = jptr;
		ljoin[dim] = jptr;
	}
	else {
		fjoin[dim] = jptr;
		ljoin[dim] = jptr;
	}
}

/**
 * @brief Adds a GridCell object to a linked list owned by Grid3D.
 * Useful for iterating over the objects and cleaning them up in the Grid3D destructor.
 * @param cptr
 * Pointer to GridCell object to add.
 */
void Grid3D::addGridCellToList(GridCell* cptr){
	if(fcell != NULL){
		lcell->next = cptr;
		lcell = cptr;
	}
	else {
		fcell = cptr;
		lcell = cptr;
	}
}

/**
 * @brief Adds a GridCell to the 3D double-linked list of GridCell objects.
 * @param newcptr
 * Pointer to the GridCell object to be added.
 */
void Grid3D::addGridCell(GridCell* newcptr){
	if(fcell != NULL){
		if(lcell->xc[0] != fcell->xc[0] + gparams->CORECELLS[0] - 1){
			newcptr->set_xcs(lcell->xc[0] + 1, lcell->xc[1], lcell->xc[2]);
			link(0, lcell, newcptr);
			link(1, locate(newcptr->xc[0], newcptr->xc[1] - 1, newcptr->xc[2]), newcptr);
			link(2, locate(newcptr->xc[0], newcptr->xc[1], newcptr->xc[2] - 1), newcptr);
		}
		else if(lcell->xc[0] == fcell->xc[0] + gparams->CORECELLS[0] - 1 && lcell->xc[1] != gparams->CORECELLS[1] - 1){
			newcptr->set_xcs(fcell->xc[0], lcell->xc[1] + 1, lcell->xc[2]);
			link(1, locate(newcptr->xc[0], newcptr->xc[1] - 1, newcptr->xc[2]), newcptr);
			link(2, locate(newcptr->xc[0], newcptr->xc[1], newcptr->xc[2] - 1), newcptr);
		}
		else if(lcell->xc[0] == fcell->xc[0] + gparams->CORECELLS[0] - 1 && lcell->xc[1] == gparams->CORECELLS[1] - 1){
			newcptr->set_xcs(fcell->xc[0], fcell->xc[1], lcell->xc[2] + 1);
			link(2, locate(newcptr->xc[0], newcptr->xc[1], newcptr->xc[2]-1), newcptr);
		}
	}
	addGridCellToList(newcptr);
}

/**
 * @brief Provides the next GridCell object along a path that traverses a 2D plane which includes the passed GridCell object.
 * @param plane
 * The 2D plane normal.
 * @param cptr
 * Pointer to previous GridCell object
 * @return Pointer to next GridCell object. If there are no more GridCell objects on this path then <code> NULL </code> is returned.
 */
GridCell* Grid3D::nextCell2D(const int& plane, GridCell* cptr){
	GridCell* newcptr = NULL;
	if((cptr->xc[(plane+2)%3]+1)%2 != 0)
		newcptr = traverse1D((plane+1)%3, 1, cptr);
	else
		newcptr = traverse1D((plane+1)%3, -1, cptr);
	if(newcptr == NULL)
		newcptr = traverse1D((plane+2)%3, 1, cptr);
	return newcptr;
}

/**
 * @brief Provides the next GridCell object along a path that traverses the entire grid.
 * @param cptr
 * Pointer to previous GridCell object
 * @return Pointer to next GridCell object. If there are no more GridCell objects on this path then <code> NULL </code> is returned.
 */
GridCell* Grid3D::nextCell3D(GridCell* cptr){
	GridCell* newcptr = nextCell2D(2, cptr);
	if(newcptr == NULL)
		newcptr = locate(fcell->xc[0], fcell->xc[1], cptr->xc[2]+1);
	return newcptr;
}

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
 * @param srcptr
 * The GridCell object that traversal in all directions starts from.
 * @param cptr
 * Pointer to previous GridCell object.
 * @param dxc
 * GridCell objects to traverse in the x direction at every step.
 * @param dyc
 * GridCell objects to traverse in the y direction at every step.
 * @param dyz
 * GridCell objects to traverse in the z direction at every step.
 * @return Pointer to next GridCell object. If there are no more GridCell objects on this path then <code> NULL </code> is returned.
 */
GridCell* Grid3D::nextSnake(GridCell* cptr, GridCell* srcptr, const int& dxc, const int& dyc, const int& dyz){
	GridCell* newcptr = traverse1D(0, dxc, cptr);
	if(newcptr == NULL && gparams->ND > 1) {
		newcptr = traverse1D(1, cptr->xc[1] - srcptr->xc[1] + dyc, srcptr);
		newcptr = traverse1D(2, cptr->xc[2] - srcptr->xc[2], newcptr);
	}
	if(newcptr == NULL && gparams->ND > 2)
		newcptr = traverse1D(2, cptr->xc[2] - srcptr->xc[2] + dyz, srcptr);
	return newcptr;
}

/**
 * @brief Provides next GridCell object along a causal path from a specified GridCell object.
 * @param cptr
 * Pointer to previous GridCell object.
 * @param srcptr
 * Pointer to the first GridCell along the path.
 * @return Pointer to next GridCell object. If there are no more GridCell objects on this path then <code> NULL </code> is returned.
 */
GridCell* Grid3D::nextCausal(GridCell* cptr, GridCell* srcptr){
	int dir[3] = {0,0,0};
	for(int i = 0; i < 3; ++i){
		if(cptr->xc[i] != srcptr->xc[i])
			dir[i] = (cptr->xc[i]-srcptr->xc[i])/std::fabs(cptr->xc[i]-srcptr->xc[i]);
		else
			dir[i] = 1;
	}
	GridCell* beginptr = srcptr;
	if(dir[0] < 0 && dir[1] > 0 && dir[2] > 0)
		beginptr = srcptr->left[0];
	if(dir[0] > 0 && dir[1] < 0 && dir[2] > 0)
		beginptr = srcptr->left[1];
	if(dir[0] > 0 && dir[1] > 0 && dir[2] < 0)
		beginptr = srcptr->left[2];
	if(dir[0] < 0 && dir[1] < 0 && dir[2] > 0)
		beginptr = srcptr->left[0]->left[1];
	if(dir[0] < 0 && dir[1] > 0 && dir[2] < 0)
		beginptr = srcptr->left[0]->left[2];
	if(dir[0] > 0 && dir[1] < 0 && dir[2] < 0)
		beginptr = srcptr->left[1]->left[2];
	if(dir[0] < 0 && dir[1] < 0 && dir[2] < 0)
		beginptr = srcptr->left[0]->left[1]->left[2];

	GridCell* newcptr = nextSnake(cptr, beginptr, dir[0], dir[1], dir[2]);

	if(newcptr == NULL && dir[0] == 1 && dir[1] == 1 && dir[2] == 1)
		newcptr = traverse3D(0, 1, 2, -1, 0, 0, srcptr);
	if(newcptr == NULL && dir[0] == -1 && dir[1] == 1 && dir[2] == 1 && gparams->ND > 1)
		newcptr = traverse3D(0, 1, 2, 0, -1, 0, srcptr);
	if(newcptr == NULL && dir[0] == 1 && dir[1] == -1 && dir[2] == 1 && gparams->ND > 1)
		newcptr = traverse3D(0, 1, 2, -1, -1, 0, srcptr);
	if(newcptr == NULL && dir[0] == -1 && dir[1] == -1 && dir[2] == 1 && gparams->ND > 2)
		newcptr = traverse3D(0, 1, 2, 0, 0, -1, srcptr);
	if(newcptr == NULL && dir[0] == 1 && dir[1] == 1 && dir[2] == -1 && gparams->ND > 2)
		newcptr = traverse3D(0, 1, 2, -1, 0, -1, srcptr);
	if(newcptr == NULL && dir[0] == -1 && dir[1] == 1 && dir[2] == -1 && gparams->ND > 2)
		newcptr = traverse3D(0, 1, 2, 0, -1, -1, srcptr);
	if(newcptr == NULL && dir[0] == 1 && dir[1] == -1 && dir[2] == -1 && gparams->ND > 2)
		newcptr = traverse3D(0, 1, 2, -1, -1, -1, srcptr);
	return newcptr;
}
void Grid3D::buildCausal(GridCell* fcausal) {
	GridCell* oldcptr = NULL;
	for (GridCell* cptr = fcausal; cptr != NULL; cptr = nextCausal(cptr, fcausal)) {
		if (oldcptr != NULL)
			oldcptr->nextcausal = cptr;
		oldcptr = cptr;
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
double Grid3D::area_join(double xj[], const int& dim) {
	double area = 0, r1, r2;
	if (gparams->GEOMETRY == CARTESIAN) {
		area = 1.0;
		for (int i = 0; i < gparams->ND; ++i) {
			if(i != dim)
				area *= dx[i];
		}
	}
	else if (gparams->GEOMETRY == CYLINDRICAL){
		/* Cylindrical polars (ND <= 2) */
		if(dim == 0) {
			r1 = (xj[0]+0.5)*dx[0];
			area = 2.0*PI*r1;
			if (gparams->ND == 2)
				area *= dx[1];
		}
		else if (dim == 1) {
			r1 = (xj[0]+0.5)*dx[0] - 0.5*dx[0];
			r2 = (xj[0]+0.5)*dx[0] + 0.5*dx[0];
			area = PI*(r2 - r1)*(r2 + r1);
		}
	}
	else if (gparams->GEOMETRY == SPHERICAL) {
		/* Spherical polars (ND = 1 only) */
		r1 = (xj[0]+0.5)*dx[0];
		area = 4.0*PI*r1*r1;
	}
	return area;
}

/**
 * @brief Calculates the volume of a GridCell given its location and the geometry of Grid3D.
 * @param cptr
 * The GridCell object for which the volume will be calculated.
 */
void Grid3D::vol_cell(GridCell* cptr) {
	double volume = 0, r1, r2;
	if(gparams->GEOMETRY == CARTESIAN){
		/* Cartesian */
		volume = 1.0;
		for (int i = 0; i < gparams->ND; ++i)
			volume *= dx[i];
	}
	else if(gparams->GEOMETRY == CYLINDRICAL){
		/* Cylindrical polars (ND <= 2) */
		r2 = (cptr->xc[0]+0.5)*dx[0] + 0.5*dx[0];
		r1 = (cptr->xc[0]+0.5)*dx[0] - 0.5*dx[0];
		volume = PI*(r2 - r1)*(r2 + r1);
		if (gparams->ND == 2)
			volume *= dx[1];
	}
	else if(gparams->GEOMETRY == SPHERICAL){
		/* Spherical polars (ND = 1 only) */
		r2 = (cptr->xc[0]+0.5)*dx[0] + 0.5*dx[0];
		r1 = (cptr->xc[0]+0.5)*dx[0] - 0.5*dx[0];
		volume = 4.0*PI*(r2 - r1)*(r2*r2 + r1*r2 + r1*r1)/3.0;
	}
	cptr->vol = volume;
}


