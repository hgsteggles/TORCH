#include "boundary.hpp"

Boundary::Boundary(int face, Grid3D* gptr) {
	this->face = face;
	this->nghosts = gptr->ORDER_S + 1;

	int dim = face%3;
	int xlength=0, ylength=0;
	if (dim == 0) {
		xlength = gptr->NCELLS[1];
		ylength = gptr->NCELLS[2];
	}
	if (dim == 1) {
		xlength = gptr->NCELLS[0];
		ylength = gptr->NCELLS[2];
	}
	if (dim == 2) {
		xlength = gptr->NCELLS[0];
		ylength = gptr->NCELLS[1];
	}
	for (int i = 0; i < xlength; i++) {
		ghostcells.push_back(std::vector<GridCell*>());
		for (int j = 0; j < ylength; j++) {
			ghostcells[i].push_back(new GridCell());
			GridCell* newghost = ghostcells[i][j];
			for (int ig = 0; ig < nghosts-1; ig++) {
				GridCell* oldghost = newghost;
				newghost = new GridCell();
				if (face < 3) {
					newghost->right[dim] = oldghost;
					oldghost->left[dim] = newghost;
				}
				else {
					newghost->left[dim] = oldghost;
					oldghost->right[dim] = newghost;
				}
			}
		}
	}
	isPartition = false;
}

Boundary::~Boundary() {
	for (int i = 0; i < (int)ghostcells.size(); i++) {
		for (int j = 0; j < (int)ghostcells[i].size(); j++) {
			GridCell* nextcptr = ghostcells[i][j];
			GridCell* cptr = nextcptr;
			while(nextcptr != NULL) {
				cptr = nextcptr;
				if (face < 3)
					nextcptr = cptr->left[face%3];
				else
					nextcptr = cptr->right[face%3];
				delete cptr;
			}
			ghostcells[i][j] = NULL;
		}
	}
}
