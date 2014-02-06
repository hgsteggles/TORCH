#include "grid3d.hpp"
#include "parameters.hpp"
#include "mpihandler.hpp"
#include "gridcell.hpp"
#include "external.hpp"
#include "partition.hpp"

Grid3D::Grid3D(const GridParameters& gp, MPIHandler& mpih) : fcell(NULL), lcell(NULL), currentTime(0.0), deltatime(0.0) {
	ND = gp.ND;
	NCELLS[0] = gp.NCELLS[0];
	NCELLS[1] = gp.NCELLS[1];
	NCELLS[2] = gp.NCELLS[2];
	TOTNCELLS[0] = gp.NCELLS[0];
	TOTNCELLS[1] = gp.NCELLS[1];
	TOTNCELLS[2] = gp.NCELLS[2];
	if (ND < 3)
		NCELLS[2] = 1;
	if (ND < 2)
		NCELLS[1] = 1;
	if (mpih.nProcessors() > 1) {
		if (NCELLS[0]%mpih.nProcessors() == 0)
			NCELLS[0] = NCELLS[0] / mpih.nProcessors();
		else {
			std::cerr << "ERROR: No. of GridCells in x dir not divisible by no. of processors." << '\n';
			exit(EXIT_FAILURE);
		}
	}
	GEOMETRY = gp.GEOMETRY;
	ORDER_S = gp.ORDER_S;
	ORDER_T = gp.ORDER_T;
	for(int i = 0; i < 3; ++i){
		fjoin[i] = NULL;
		ljoin[i] = NULL;
	}
	dx[0] = 1.0/(double)gp.NCELLS[0];
	dx[1] = 1.0/(double)gp.NCELLS[1];
	dx[2] = 1.0/(double)gp.NCELLS[2];

	for(int i = 0, prc = 0; i < (NCELLS[0]*NCELLS[1]*NCELLS[2]); ++i){
		int prcnow = (int) (100*i/(double)(NCELLS[0]*NCELLS[1]*NCELLS[2]));
		if(prcnow - prc >= 10) {
			std::cout << "GRID3D: Adding GridCell objects... " << prcnow << "%" << '\n';
			prc = prcnow;
		}
		GridCell* newcptr = new GridCell();
		if (i == 0) {
			newcptr->set_xcs(mpih.getRank()*NCELLS[0], 0, 0);
		}
		addGridCell(newcptr);
		vol_cell(newcptr);
	}
	std::cout << mpih.cname() <<  "GRID3D: set up with " << NCELLS[0]*NCELLS[1]*NCELLS[2] << " cells." << '\n';

	if(ND > 0) {
		if (mpih.getRank() == 0) {
			boundaries.push_back(new ExternalBoundary(0, REFLECTING, this));
			std::cout << mpih.cname() << "EXTBOUNDARY: set up with " << boundaries[0]->nghosts << " ghostcells on face " << 0 << '\n';
		}
		else {
			int destination = mpih.getRank()-1;
			boundaries.push_back(new Partition(0, this, destination, mpih));
			std::cout << mpih.cname() << "PARTITION: set up with " << boundaries[0]->nghosts << " ghostcells on face " << 0 << '\n';
		}
		boundaryLink(boundaries[0]);
		if (mpih.getRank() == mpih.nProcessors()-1) {
			boundaries.push_back(new ExternalBoundary(3, FREE, this));
			std::cout << mpih.cname() << "EXTBOUNDARY: set up with " << boundaries[1]->nghosts << " ghostcells on face " << 3 << '\n';
		}
		else {
			boundaries.push_back(new Partition(3, this, mpih.getRank()+1, mpih));
			std::cout << mpih.cname() << "PARTITION: set up with " << boundaries[1]->nghosts << " ghostcells on face " << 3 << '\n';
		}
		boundaryLink(boundaries[1]);
	}
	if(ND > 1) {
		boundaries.push_back(new ExternalBoundary(1, REFLECTING, this));
		boundaryLink(boundaries[2]);
		std::cout << mpih.cname() << "EXTBOUNDARY: set up with " << boundaries[2]->nghosts << " ghostcells on face " << 1 << '\n';
		boundaries.push_back(new ExternalBoundary(4, FREE, this));
		boundaryLink(boundaries[3]);
		std::cout << mpih.cname() << "EXTBOUNDARY: set up with " << boundaries[3]->nghosts << " ghostcells on face " << 4 << '\n';
	}
	if(ND > 2) {
		boundaries.push_back(new ExternalBoundary(2, REFLECTING, this));
		boundaryLink(boundaries[4]);
		std::cout << mpih.cname() << "EXTBOUNDARY: set up with " << boundaries[4]->nghosts << " ghostcells on face " << 2 << '\n';
		boundaries.push_back(new ExternalBoundary(5, FREE, this));
		boundaryLink(boundaries[5]);
		std::cout << mpih.cname() << "EXTBOUNDARY: set up with " << boundaries[5]->nghosts << " ghostcells on face " << 5 << '\n';
	}
	std::cout << "printinfo 1\n";
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
	for(int dim = 0; dim < ND; ++dim){
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
			for (int i = 0; i < 3; ++i)
				bcptr->xc[i] = cptr->xc[i];
			bcptr->xc[dim] = cptr->xc[dim] - 1;
			weakLink(dim, bcptr, cptr);
		}
	}
	else {
		for (GridCell* cptr = traverse1D(dim, NCELLS[dim]-1, fcell); cptr != NULL; cptr = nextCell2D(dim, cptr)) {
			GridCell* bcptr = bptr->ghostcells[cptr->xc[ii]][cptr->xc[jj]];
			for (int i = 0; i < 3; ++i)
				bcptr->xc[i] = cptr->xc[i];
			bcptr->xc[dim] = cptr->xc[dim] + 1;
			weakLink(dim, cptr, bcptr);
		}
	}
}
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
GridCell* Grid3D::traverse3D(const int& d1, const int& d2, const int& d3, const int& dc1, const int& dc2, const int& dc3, GridCell* cptr){
	return traverse1D(d3, dc3, traverse1D(d2, dc2, traverse1D(d1, dc1, cptr)));
}
GridCell* Grid3D::traverseOverJoins1D(const int& dim, const int& dxc, GridCell* cptr){
	GridCell* newcptr = cptr;
	for(int i = 0; i < abs(dxc) && newcptr != NULL; ++i){
		if (dxc > 0 && newcptr->rjoin[dim] != NULL)
			newcptr = newcptr->rjoin[dim]->rcell;
		else if (dxc <= 0 && newcptr->ljoin[dim] != NULL)
			newcptr = newcptr->ljoin[dim]->lcell;
		else
			newcptr = NULL;
	}
	return newcptr;
}
GridCell* Grid3D::traverseOverJoins3D(const int& d1, const int& d2, const int& d3, const int& dc1, const int& dc2, const int& dc3, GridCell* cptr){
	return traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d1, dc1, cptr)));
}
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
void Grid3D::addGridCell(GridCell* newcptr){
	if(fcell != NULL){
		if(lcell->xc[0] != NCELLS[0] - 1){
			newcptr->set_xcs(lcell->xc[0] + 1, lcell->xc[1], lcell->xc[2]);
			link(0, lcell, newcptr);
			link(1, locate(newcptr->xc[0], newcptr->xc[1] - 1, newcptr->xc[2]), newcptr);
			link(2, locate(newcptr->xc[0], newcptr->xc[1], newcptr->xc[2] - 1), newcptr);
		}
		else if(lcell->xc[0] == NCELLS[0] - 1 && lcell->xc[1] != NCELLS[1] - 1){
			newcptr->set_xcs(fcell->xc[0], lcell->xc[1] + 1, lcell->xc[2]);
			link(1, locate(newcptr->xc[0], newcptr->xc[1] - 1, newcptr->xc[2]), newcptr);
			link(2, locate(newcptr->xc[0], newcptr->xc[1], newcptr->xc[2] - 1), newcptr);
		}
		else if(lcell->xc[0] == NCELLS[0] - 1 && lcell->xc[1] == NCELLS[1] - 1){
			newcptr->set_xcs(fcell->xc[0], fcell->xc[1], lcell->xc[2] + 1);
			link(2, locate(newcptr->xc[0], newcptr->xc[1], newcptr->xc[2]-1), newcptr);
		}
	}
	addGridCellToList(newcptr);
}

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

GridCell* Grid3D::nextCell3D(GridCell* cptr){
	GridCell* newcptr = nextCell2D(2, cptr);
	if(newcptr == NULL)
		newcptr = locate(fcell->xc[0], fcell->xc[1], cptr->xc[2]+1);
	return newcptr;
}
GridCell* Grid3D::nextSnake(GridCell* cptr, GridCell* srcptr, const int& dxc, const int& dyc, const int& dyz){
	GridCell* newcptr = traverse1D(0, dxc, cptr);
	if(newcptr == NULL && ND > 1) {
		newcptr = traverse1D(1, cptr->xc[1] - srcptr->xc[1] + dyc, srcptr);
		newcptr = traverse1D(2, cptr->xc[2] - srcptr->xc[2], newcptr);
	}
	if(newcptr == NULL && ND > 2)
		newcptr = traverse1D(2, cptr->xc[2] - srcptr->xc[2] + dyz, srcptr);
	return newcptr;
}
GridCell* Grid3D::nextCausal(GridCell* cptr, GridCell* srcptr){
	int dir[3] = {0,0,0};
	for(int i = 0; i < 3; ++i){
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

	GridCell* newcptr = nextSnake(cptr, beginptr, dir[0], dir[1], dir[2]);
	if (cptr->xc[0] == 127 && cptr->xc[1] == 0 && cptr->xc[2] == 1) {
		if (newcptr != NULL)
			std::cout << newcptr->xc[0] << " " << newcptr->xc[1] << " " << newcptr->xc[2] << '\n';
	}
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
void Grid3D::buildCausal(GridCell* fcausal) {
	GridCell* oldcptr = NULL;
	for (GridCell* cptr = fcausal; cptr != NULL; cptr = nextCausal(cptr, fcausal)) {
		if (oldcptr != NULL)
			oldcptr->nextcausal = cptr;
		oldcptr = cptr;
	}
}
double Grid3D::area_join(double xj[], const int& dim) {
	double area = 0, r1, r2;
	if (GEOMETRY == CARTESIAN) {
		area = 1.0;
		for (int i = 0; i < ND; ++i) {
			if(i != dim)
				area *= dx[i];
		}
	}
	else if (GEOMETRY == CYLINDRICAL){
		/* Cylindrical polars (ND <= 2) */
		if(dim == 0) {
			r1 = (xj[0]+0.5)*dx[0];
			area = 2.0*PI*r1;
			if (ND == 2)
				area *= dx[1];
		}
		else if (dim == 1) {
			r1 = (xj[0]+0.5)*dx[0] - 0.5*dx[0];
			r2 = (xj[0]+0.5)*dx[0] + 0.5*dx[0];
			area = PI*(r2 - r1)*(r2 + r1);
		}
	}
	else if (GEOMETRY == SPHERICAL) {
		/* Spherical polars (ND = 1 only) */
		r1 = (xj[0]+0.5)*dx[0];
		area = 4.0*PI*r1*r1;
	}
	return area;
}
void Grid3D::vol_cell(GridCell* cptr) {
	double volume = 0, r1, r2;
	if(GEOMETRY == CARTESIAN){
		/* Cartesian */
		volume = 1.0;
		for (int i = 0; i < ND; ++i)
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


