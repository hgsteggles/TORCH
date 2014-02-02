/*
 * partition.cpp
 *
 *  Created on: 17 Jan 2014
 *      Author: "Harrison Steggles"
 */

#include "partition.hpp"

Partition::Partition(int face, Grid3D* gptr, int dest,  MPIHandler& mpih) : Boundary(face, gptr), destination(dest), mpihandler(mpih) {
	isPartition = true;
}

void Partition::applyBC() {
	int dim = face%3;
	unsigned int NELEMENTS = ghostcells.size()*ghostcells[0].size()*nghosts*NU;
	double msgArray[NELEMENTS];
	int id = 0;
	for (int i = 0; i < (int)ghostcells.size(); i++) {
		for (int j = 0; j < (int)ghostcells[i].size(); j++) {
			GridCell* cptr = NULL;
			if (face < 3)
				cptr = ghostcells[i][j]->rjoin[dim]->rcell;
			else
				cptr = ghostcells[i][j]->ljoin[dim]->lcell;
			for (int k = 0; k < nghosts && cptr != NULL; k++) {
				for(int iu = 0; iu < NU; iu++) {
					msgArray[id] = cptr->Q[iu];
					id++;
				}
				if(face < 3)
					cptr = cptr->rjoin[dim]->rcell;
				else
					cptr = cptr->ljoin[dim]->lcell;
			}
		}
	}
	mpihandler.send(destination, PARTITION_MSG, msgArray, NELEMENTS);
	mpihandler.receive(destination, PARTITION_MSG, msgArray, NELEMENTS);
	mpihandler.wait();
	id = 0;
	for (int i = 0; i < (int)ghostcells.size(); i++) {
		for (int j = 0; j < (int)ghostcells[i].size(); j++) {
			GridCell* ghost = ghostcells[i][j];
			for (int k = 0; k < nghosts; k++) {
				for(int iu = 0; iu < NU; iu++) {
					ghost->Q[iu] = msgArray[id];
					id++;
				}
				if(face < 3)
					ghost = ghost->left[dim];
				else
					ghost = ghost->right[dim];
			}
		}
	}

}
/*
void Partition::applyBC() {
	int dim = face%3;
	int tag = 0;
	for (int i = 0; i < (int)ghostcells.size(); i++) {
		for (int j = 0; j < (int)ghostcells[i].size(); j++) {
			GridCell* ghost = ghostcells[i][j];
			GridCell* cptr = NULL;
			if (face < 3)
				cptr = ghostcells[i][j]->rjoin[dim]->rcell;
			else
				cptr = ghostcells[i][j]->ljoin[dim]->lcell;

			for (int k = 0; k < nghosts && cptr != NULL; k++) {
				for(int iu = 0; iu < NU; iu++) {
					mpihandler.send(destination, tag, cptr->Q[iu]);
					tag++;
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
	tag = 0;
	for (int i = 0; i < (int)ghostcells.size(); i++) {
		for (int j = 0; j < (int)ghostcells[i].size(); j++) {
			GridCell* ghost = ghostcells[i][j];
			while (ghost != NULL) {
				for(int iu = 0; iu < NU; iu++) {
					mpihandler.receive(destination, tag, ghost->Q[iu]);
					tag++;
				}
				if(face < 3)
					ghost = ghost->left[dim];
				else
					ghost = ghost->right[dim];
			}
		}
	}
	mpihandler.wait();
}
*/
