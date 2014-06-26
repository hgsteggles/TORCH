/**
 * @file boundary.cpp
 */

#include "boundary.h"
#include "grid3d.h"
#include "gridcell.h"
#include "parameters.h"

#include <stddef.h>

/**
 * @brief Boundary constructor.
 * Requires a face and pointer to Grid3D so that the size of the 2D GridCell vector to be built is known.
 * @param face An integer specifying the face this Boundary is attached to (face%3 gives dimension, left
 * boundaries have face<3 and right boundaries have face>=3).
 * @param gptr Pointer to Grid3D instance that this Boundary is to be attached to.
 */
Boundary::Boundary(const int& face, const int nOfGhosts, Grid3D* gptr) {
	this->face = face;
	this->nghosts = nOfGhosts;

	int dim = face%3;
	int xlength = 0, ylength = 0;
	if (dim == 0) {
		xlength = gptr->gparams->CORECELLS[1];
		ylength = gptr->gparams->CORECELLS[2];
	}
	if (dim == 1) {
		xlength = gptr->gparams->CORECELLS[0];
		ylength = gptr->gparams->CORECELLS[2];
	}
	if (dim == 2) {
		xlength = gptr->gparams->CORECELLS[0];
		ylength = gptr->gparams->CORECELLS[1];
	}
	for (int i = 0; i < xlength; ++i) {
		ghostcells.push_back(std::vector<GridCell*>());
		for (int j = 0; j < ylength; ++j) {
			ghostcells[i].push_back(new GridCell());
			GridCell* newghost = ghostcells[i][j];
			for (int ig = 0; ig < nghosts-1; ++ig) {
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

/**
 * @brief Boundary destructor.
 * Deletes all GridCell objects allocated
 */
Boundary::~Boundary() {
	for (int i = 0; i < (int)ghostcells.size(); ++i) {
		for (int j = 0; j < (int)ghostcells[i].size(); ++j) {
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
