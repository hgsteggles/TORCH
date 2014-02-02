#include "external.hpp"

ExternalBoundary::ExternalBoundary(int face, Condition bcond, Grid3D* gptr) : Boundary(face, gptr), bc(bcond) { }

void ExternalBoundary::applyBC() {
	int dim = face%3;
	for (int i = 0; i < (int)ghostcells.size(); i++) {
		for (int j = 0; j < (int)ghostcells[i].size(); j++) {
			GridCell* ghost = ghostcells[i][j];
			GridCell* cptr = NULL;
			if (face < 3)
				cptr = ghostcells[i][j]->rjoin[dim]->rcell;
			else
				cptr = ghostcells[i][j]->ljoin[dim]->lcell;
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
						ghost->Q[j] = cptr->Q[j];
					if((cptr->Q[ivel+dim] < 0 && face < 3) || (cptr->Q[ivel+dim] > 0 && face >= 3))
						ghost->Q[ivel+dim] = -1.0*cptr->Q[ivel+dim];
				}
				if(bc == FREE){
					for(int iu = 0; iu < NU; iu++)
						ghost->Q[j] = cptr->Q[j];
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
}


