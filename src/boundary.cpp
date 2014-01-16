/*
 * boundary.cpp
 *
 *  Created on: 14 Jan 2014
 *      Author: harry
 */

#include "boundary.hpp"

Boundary::Boundary(int face, Condition bc, Grid3D* gptr) {
	this->bc = bc;
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

void Boundary::applyBC(){
	int dim = face%3;
	for (int i = 0; i < (int)ghostcells.size(); i++) {
		for (int j = 0; j < (int)ghostcells[i].size(); j++) {
			GridCell* ghost = ghostcells[i][j];
			GridCell* cptr = NULL;
			if (face < 3)
				cptr = ghostcells[i][j]->rjoin[face%3]->rcell;
			else
				cptr = ghostcells[i][j]->ljoin[face%3]->lcell;
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
					cptr = cptr->rjoin[face%3]->rcell;
					ghost = ghost->left[face%3];
				}
				else {
					cptr = cptr->ljoin[face%3]->lcell;
					ghost = ghost->right[face%3];
				}
			}
		}
	}
}
