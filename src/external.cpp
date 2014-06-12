/**
 * @file external.cpp
 */

#include "external.hpp"
#include "grid3d.hpp"
#include "gridcell.hpp"


#include <stddef.h> // NULL
#include <stdlib.h> // exit
#include <iostream>

/**
 * @brief ExternalBoundary constructor.
 * @param face
 * The face of the grid that the boundary is connected to (face < 3: left face, else: right face).
 * @param bcond
 * The boundary condition that must be applied to this boundary at the start of each hydro/radiation update step.
 * @param gptr
 * A pointer to the Grid3D object that the Boundary should be linked to with Grid3D::boundaryLink(Boundary*).
 */
ExternalBoundary::ExternalBoundary(const int face, const int nOfGhosts, const Condition& bcond, Grid3D* gptr) :
	Boundary(face, nOfGhosts, gptr),
	bc(bcond) {
}

/**
 * @brief Applies the boundary condition that was passed to this object through its constructor.
 */
void ExternalBoundary::applyBC() {
	int dim = face%3;
	for (int i = 0; i < (int)ghostcells.size(); ++i) {
		for (int j = 0; j < (int)ghostcells[i].size(); ++j) {
			GridCell* ghost = ghostcells[i][j];
			GridCell* cptr = NULL;
			if (face < 3)
				cptr = ghostcells[i][j]->rjoin[dim]->rcell;
			else{
				if (ghostcells[i][j]->ljoin[dim] == NULL) {
					ghostcells[i][j]->printInfo();
					exit(EXIT_FAILURE);
				}
				else
					cptr = ghostcells[i][j]->ljoin[dim]->lcell;
			}
			while (ghost != NULL && cptr != NULL) {
				if(bc == REFLECTING){
					for(int iu = 0; iu < NU; iu++)
						ghost->Q[iu] = cptr->Q[iu];
					ghost->Q[ivel+dim] = -cptr->Q[ivel+dim];
				}
				if(bc == OUTFLOW){
					for(int iu = 0; iu < NU; iu++)
						ghost->Q[iu] = cptr->Q[iu];
					if((cptr->Q[ivel+dim] > 0 && face < 3) || (cptr->Q[ivel+dim] < 0 && face >= 3))
						ghost->Q[ivel+dim] = -1.0*cptr->Q[ivel+dim];
				}
				if(bc == INFLOW){
					for(int iu = 0; iu < NU; iu++)
						ghost->Q[iu] = cptr->Q[iu];
					if((cptr->Q[ivel+dim] < 0 && face < 3) || (cptr->Q[ivel+dim] > 0 && face >= 3))
						ghost->Q[ivel+dim] = -1.0*cptr->Q[ivel+dim];
				}
				if(bc == FREE){
					for(int iu = 0; iu < NU; iu++)
						ghost->Q[iu] = cptr->Q[iu];
				}
				if(face < 3) {
					cptr = cptr->rjoin[dim]->rcell;
					ghost = ghost->left[dim];
				}
				else {
					cptr = cptr->ljoin[dim]->lcell;
					ghost = ghost->right[dim];
				}
			}
		}
	}
	for (int i = 0; i < (int)ghostcells.size(); ++i) {
		for (int j = 0; j < (int)ghostcells[i].size(); ++j) {
			GridCell* ghost = ghostcells[i][j];
			GridCell* cptr = NULL;
			if (face < 3)
				cptr = ghostcells[i][j]->rjoin[dim]->rcell;
			else
				cptr = ghostcells[i][j]->ljoin[dim]->lcell;
			for (int ir = 0; ir < NR; ++ir)
				ghost->R[ir] = cptr->R[ir];
		}
	}
}


