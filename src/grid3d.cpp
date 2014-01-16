#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "grid3d.hpp"

using namespace std;

Grid3D::Grid3D(const GridParameters& gp) : fcell(NULL), lcell(NULL), fsrc(NULL), fboundary(NULL) {
	ND = gp.ND;
	NCELLS[0] = gp.NCELLS[0];
	NCELLS[1] = gp.NCELLS[1];
	NCELLS[2] = gp.NCELLS[2];
	GEOMETRY = gp.GEOMETRY;
	ORDER_S = gp.ORDER_S;
	ORDER_T = gp.ORDER_T;
	GridCell *newcptr;
	for(int i = 0; i < 3; i++)
		fjoin[i] = NULL;
	dx[0] = 1.0/(double)NCELLS[0];
	dx[1] = 1.0/(double)NCELLS[1];
	dx[2] = 1.0/(double)NCELLS[2];


	for(int i=0; i < (NCELLS[0]*NCELLS[1]*NCELLS[2]); i++){
		newcptr = new GridCell();
		addGridCell(newcptr);
		vol_cell(newcptr);
	}

	cout << "MAIN: grid set up with " << NCELLS[0]*NCELLS[1]*NCELLS[2] << " cells." << '\n';

	if(ND > 0) {
		boundaries.push_back(new Boundary(0, REFLECTING, this));
		boundaryLink(boundaries[0]);
		cout << "MAIN: Boundary set up with " << boundaries[0]->nghosts << " ghostcells." << endl;
		boundaries.push_back(new Boundary(3, REFLECTING, this));
		boundaryLink(boundaries[1]);
		cout << "MAIN: Boundary set up with " << boundaries[1]->nghosts << " ghostcells." << endl;
	}
	if(ND > 1) {
		boundaries.push_back(new Boundary(1, FREE, this));
		boundaryLink(boundaries[2]);
		cout << "MAIN: Boundary set up with " << boundaries[2]->nghosts << " ghostcells." << endl;
		boundaries.push_back(new Boundary(4, FREE, this));
		boundaryLink(boundaries[3]);
		cout << "MAIN: Boundary set up with " << boundaries[3]->nghosts << " ghostcells." << endl;
	}
	if(ND > 2) {
		boundaries.push_back(new Boundary(2, FREE, this));
		boundaryLink(boundaries[4]);
		cout << "MAIN: Boundary set up with " << boundaries[4]->nghosts << " ghostcells." << endl;
		boundaries.push_back(new Boundary(5, FREE, this));
		boundaryLink(boundaries[5]);
		cout << "MAIN: Boundary set up with " << boundaries[5]->nghosts << " ghostcells." << endl;
	}

	s_total++;
}
Grid3D::~Grid3D(){ 
	GridJoin* nextjptr = NULL;
	GridCell* nextcptr = NULL;
	for (int i = 0; i < (int)boundaries.size(); i++) {
		delete boundaries[i];
		boundaries[i] = NULL;
	}
	for(int dim = 0; dim < 3; dim++){
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

void Grid3D::dice(){
	for(GridCell* cptr = fcell; cptr != NULL; cptr = cptr->next){
		for(int dim = 0; dim < ND; dim++){
			if(cptr->xc[dim]%2 == 0)
				cptr->right[dim] = NULL;
		}
	}
}

GridCell* Grid3D::start() {
	return fcell;
}
void Grid3D::link(int dim, GridCell* lcptr, GridCell* rcptr){
	if(lcptr != NULL && rcptr != NULL){
		rcptr->left[dim] = lcptr;
		lcptr->right[dim] = rcptr;
		GridJoin* jptr = new GridJoin();
		rcptr->ljoin[dim] = jptr;
		lcptr->rjoin[dim] = jptr;
		jptr->rcell = rcptr;
		jptr->lcell = lcptr;
		for(int i = 0; i < 3; i++)
			jptr->xj[i] = jptr->rcell->xc[i];
		jptr->xj[dim] = jptr->rcell->xc[dim] - 0.5;
		jptr->area = area_join(jptr->xj, dim);
		addJoinToList(jptr, dim);
	}
}
void Grid3D::weakLink(int dim, GridCell* lcptr, GridCell* rcptr){
	if(lcptr != NULL && rcptr != NULL){
		GridJoin* jptr = new GridJoin();
		rcptr->ljoin[dim] = jptr;
		lcptr->rjoin[dim] = jptr;
		jptr->rcell = rcptr;
		jptr->lcell = lcptr;
		for(int i = 0; i < 3; i++)
			jptr->xj[i] = jptr->rcell->xc[i];
		jptr->xj[dim] = jptr->rcell->xc[dim] - 0.5;
		jptr->area = area_join(jptr->xj, dim);
		addJoinToList(jptr, dim);
	}
}
void Grid3D::boundaryLink(Boundary* bptr) {
	int dim = bptr->face%3;
	int ii = 0;
	int jj = 2;
	if (dim == 0)
		ii = 1;
	if (dim == 2)
		jj = 1;

	if (bptr->face < 3) {
		for (GridCell* cptr = fcell; cptr != NULL; cptr = nextCell2D(dim, cptr)) {
			GridCell* bcptr = bptr->ghostcells[cptr->xc[ii]][cptr->xc[jj]];
			for (int i = 0; i < 3; i++)
				bcptr->xc[i] = cptr->xc[i];
			bcptr->xc[dim] = cptr->xc[dim] - 1;
			weakLink(dim, bcptr, cptr);
		}
	}
	else {
		for (GridCell* cptr = traverse1D(dim, NCELLS[dim]-1, fcell); cptr != NULL; cptr = nextCell2D(dim, cptr)) {
			GridCell* bcptr = bptr->ghostcells[cptr->xc[ii]][cptr->xc[jj]];
			for (int i = 0; i < 3; i++)
				bcptr->xc[i] = cptr->xc[i];
			bcptr->xc[dim] = cptr->xc[dim] + 1;
			weakLink(dim, cptr, bcptr);
		}
	}
}
GridCell* Grid3D::traverse1D(int dim, int dxc, GridCell* cptr){
	GridCell* newcptr = cptr;
	for(int i = 0; i < abs(dxc) && newcptr != NULL; i++){
		if(dxc > 0)
			newcptr = newcptr->right[dim];
		else
			newcptr = newcptr->left[dim]; 
	}
	return newcptr;
}
GridCell* Grid3D::traverse3D(int d1, int d2, int d3, int dc1, int dc2, int dc3, GridCell* cptr){
	return traverse1D(d3, dc3, traverse1D(d2, dc2, traverse1D(d1, dc1, cptr)));
}
GridCell* Grid3D::locate(int x, int y, int z){
	GridCell* cptr = fcell;
	for(int i = 0, tc[3] = {x, y, z}; i < 3; i++){
		for(int j = 0; (cptr != NULL) && (cptr->get_xc(i) != tc[i]); j++)
			cptr = cptr->right[i];
	}
	return cptr;
}
void Grid3D::addSource(int x, int y, int z){
	GridCell* nextsrc = locate(x, y, z);
	if(fsrc != NULL){
		GridCell* lsrc = fsrc;
		while(lsrc->next != NULL)
			lsrc = lsrc->next;
		lsrc->next = nextsrc;
	}
	else
		fsrc = nextsrc;
}
void Grid3D::addJoinToList(GridJoin* jptr, int dim){
	if(fjoin[dim] != NULL){
		GridJoin* ljoin = fjoin[dim];
		while(ljoin->next != NULL)
			ljoin = ljoin->next;
		ljoin->next = jptr;
	}
	else
		fjoin[dim] = jptr;
}
void Grid3D::addGridCellToList(GridCell* cptr){
	if(fcell != NULL){
		GridCell* endcell = fcell;
		while(endcell->next != NULL)
			endcell = endcell->next;
		endcell->next = cptr;
	}
	else
		fcell = cptr;
}
void Grid3D::deleteBoundary(Boundary* bptr) {
	delete bptr;
	bptr = NULL;
}
void Grid3D::addGridCell(GridCell* newcptr){
	if(fcell != NULL){
		if(lcell->get_xc(0) != NCELLS[0] - 1){
			newcptr->set_xcs(lcell->get_xc(0) + 1, lcell->get_xc(1), lcell->get_xc(2));
			link(0, lcell, newcptr);
			link(1, locate(newcptr->get_xc(0), newcptr->get_xc(1) - 1, newcptr->get_xc(2)), newcptr);
			link(2, locate(newcptr->get_xc(0), newcptr->get_xc(1), newcptr->get_xc(2) - 1), newcptr);
		}
		else if(lcell->get_xc(0) == NCELLS[0] - 1 && lcell->get_xc(1) != NCELLS[1] - 1){
			newcptr->set_xcs(0, lcell->get_xc(1) + 1, lcell->get_xc(2));
			link(1, locate(newcptr->get_xc(0), newcptr->get_xc(1) - 1, newcptr->get_xc(2)), newcptr);
			link(2, locate(newcptr->get_xc(0), newcptr->get_xc(1), newcptr->get_xc(2) - 1), newcptr);
		}
		else if(lcell->get_xc(0) == NCELLS[0] - 1 && lcell->get_xc(1) == NCELLS[1] - 1){
			newcptr->set_xcs(0, 0, (lcell->get_xc(2) )+1);
			link(2, locate(newcptr->get_xc(0), newcptr->get_xc(1), newcptr->get_xc(2)-1), newcptr);
		}
	}
	lcell = newcptr;
	addGridCellToList(newcptr);
}

GridCell* Grid3D::nextCell2D(int plane, GridCell* cptr){
	GridCell* newcptr = NULL;
	if((cptr->xc[(plane+2)%3]+1)%2 != 0)
		newcptr = traverse1D((plane+1)%3, 1, cptr);
	else
		newcptr = traverse1D((plane+1)%3, -1, cptr);
	if(newcptr == NULL)
		newcptr = traverse1D((plane+2)%3, 1, cptr);
	return newcptr;
}

GridCell* Grid3D::nextCell3D(GridCell* cptr){
	GridCell* newcptr = nextCell2D(2, cptr);
	if(newcptr == NULL)
		newcptr = locate(0, 0, cptr->xc[2]+1);
	return newcptr;
}
GridJoin* Grid3D::nextGridlink(int dim, GridJoin* jptr){
	GridJoin* newjptr = NULL;
	if(jptr->rcell->get_xc(0) != -1){
		if(nextCell3D(jptr->rcell) != NULL)
			newjptr = nextCell3D(jptr->rcell)->ljoin[dim];
		if(newjptr == NULL)
			newjptr = traverse1D(dim, NCELLS[dim]-1, locate(0, 0, 0))->rjoin[dim];
	}
	else{
		if(nextCell2D(dim, newjptr->lcell) != NULL)
			newjptr = nextCell2D(dim, newjptr->lcell)->rjoin[dim];
	}
	return newjptr;
}
GridCell* Grid3D::nextSnake(GridCell* srcptr, GridCell* cptr, int dxc, int dyc, int dyz){
	GridCell* newcptr = traverse1D(0, dxc, cptr);
	if(newcptr == NULL && ND > 1)
		newcptr = traverse1D(1, cptr->xc[1] - srcptr->xc[1] + dyc, srcptr);
	if(newcptr == NULL && ND > 2)
		newcptr = traverse1D(2, cptr->xc[2] - srcptr->xc[2] + dyz, srcptr);
	return newcptr;
}
GridCell* Grid3D::nextCausal(GridCell* cptr, GridCell* srcptr){
	int dir[3] = {0,0,0};
	for(int i = 0; i < 3; i++){
		if(cptr->xc[i] != srcptr->xc[i])
			dir[i] = (cptr->xc[i]-srcptr->xc[i])/fabs(cptr->xc[i]-srcptr->xc[i]);
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

	GridCell* newcptr = nextSnake(beginptr, cptr, dir[0], dir[1], dir[2]);

	if(newcptr == NULL && dir[0] == 1 && dir[1] == 1 && dir[2] == 1)
		newcptr = traverse3D(0, 1, 2, -1, 0, 0, srcptr);
	if(newcptr == NULL && dir[0] == -1 && dir[1] == 1 && dir[2] == 1 && ND > 1)
		newcptr = traverse3D(0, 1, 2, 0, -1, 0, srcptr);
	if(newcptr == NULL && dir[0] == 1 && dir[1] == -1 && dir[2] == 1 && ND > 1)
		newcptr = traverse3D(0, 1, 2, -1, -1, 0, srcptr);
	if(newcptr == NULL && dir[0] == -1 && dir[1] == -1 && dir[2] == 1 && ND > 2)
		newcptr = traverse3D(0, 1, 2, 0, 0, -1, srcptr);
	if(newcptr == NULL && dir[0] == 1 && dir[1] == 1 && dir[2] == -1 && ND > 2)
		newcptr = traverse3D(0, 1, 2, -1, 0, -1, srcptr);
	if(newcptr == NULL && dir[0] == -1 && dir[1] == 1 && dir[2] == -1 && ND > 2)
		newcptr = traverse3D(0, 1, 2, 0, -1, -1, srcptr);
	if(newcptr == NULL && dir[0] == 1 && dir[1] == -1 && dir[2] == -1 && ND > 2)
		newcptr = traverse3D(0, 1, 2, -1, -1, -1, srcptr);
	return newcptr;
}
double Grid3D::area_join(double xj[], int dim){
	double area = 0, r1, r2;
	if(GEOMETRY == CARTESIAN){
		area = 1.0;
		for (int i = 0; i < ND; i++){
			if(i != dim)
				area *= dx[i];
		}
	}
	else if (GEOMETRY == CYLINDRICAL){
		/* Cylindrical polars (ND <= 2) */
		if (dim == 0){
			r1 = (xj[0]+0.5)*dx[0];
			area = 2.0*PI*r1;
			if (ND == 2)
				area *= dx[1];
		}
		else if(dim == 1){
			r1 = (xj[0]+0.5)*dx[0] - 0.5*dx[0];
			r2 = (xj[0]+0.5)*dx[0] + 0.5*dx[0];
			area = PI*(r2 - r1)*(r2 + r1);
		}
	}
	else if (GEOMETRY == SPHERICAL){
		/* Spherical polars (ND = 1 only) */
		r1 = (xj[0]+0.5)*dx[0];
		area = 4.0*PI*r1*r1;
	}
	return area;
}
void Grid3D::vol_cell(GridCell* cptr){
	double volume = 0, r1, r2;
	if(GEOMETRY == CARTESIAN){
		/* Cartesian */
		volume = 1.0;
		for (int i = 0; i < ND; i++)
			volume *= dx[i];
	}
	else if(GEOMETRY == CYLINDRICAL){
		/* Cylindrical polars (ND <= 2) */
		r2 = (cptr->xc[0]+0.5)*dx[0] + 0.5*dx[0];
		r1 = (cptr->xc[0]+0.5)*dx[0] - 0.5*dx[0];
		volume = PI*(r2 - r1)*(r2 + r1);
		if (ND == 2)
			volume *= dx[1];
	}
	else if(GEOMETRY == SPHERICAL){
		/* Spherical polars (ND = 1 only) */
		r2 = (cptr->xc[0]+0.5)*dx[0] + 0.5*dx[0];
		r1 = (cptr->xc[0]+0.5)*dx[0] - 0.5*dx[0];
		volume = 4.0*PI*(r2 - r1)*(r2*r2 + r1*r2 + r1*r1)/3.0;
	}
	cptr->vol = volume;
}


