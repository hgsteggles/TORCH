/*
 * partition.cpp
 *
 *  Created on: 17 Jan 2014
 *      Author: "Harrison Steggles"
 */

#include "partition.hpp"
#include "boundary.hpp"
#include "mpihandler2.hpp"
#include "constants.hpp"
#include "gridcell.hpp"

#include <fstream>

Partition::Partition(const int face, Grid3D* gptr, const int& dest,  MPIHandler& mpih) : Boundary(face, gptr), destination(dest), mpihandler(mpih) {
	isPartition = true;
}

void Partition::applyBC() {
	int dim = face%3;
	const unsigned int NELEMENTS = ghostcells.size()*ghostcells[0].size()*nghosts*NU;
	double* msgArray = new double[NELEMENTS];
	int id = 0;
	for (int i = 0; i < (int)ghostcells.size(); ++i) {
		for (int j = 0; j < (int)ghostcells[i].size(); ++j) {
			GridCell* cptr = NULL;
			if (face < 3)
				cptr = ghostcells[i][j]->rjoin[dim]->rcell;
			else
				cptr = ghostcells[i][j]->ljoin[dim]->lcell;
			for (int k = 0; k < nghosts && cptr != NULL; ++k) {
				for(int iu = 0; iu < NU; ++iu) {
					msgArray[id] = cptr->Q[iu];
					++id;
				}
				if(face < 3)
					cptr = cptr->rjoin[dim]->rcell;
				else
					cptr = cptr->ljoin[dim]->lcell;
			}
		}
	}
	
	mpihandler.exchange(msgArray, NELEMENTS, destination, PARTITION_MSG);

	id = 0;
	for (int i = 0; i < (int)ghostcells.size(); ++i) {
		for (int j = 0; j < (int)ghostcells[i].size(); ++j) {
			GridCell* ghost = ghostcells[i][j];
			for (int k = 0; k < nghosts; ++k) {
				for(int iu = 0; iu < NU; ++iu) {
					ghost->Q[iu] = msgArray[id];
					++id;
				}
				if(face < 3)
					ghost = ghost->left[dim];
				else
					ghost = ghost->right[dim];
			}
		}
	}
	delete[] msgArray;
}
/*
void Partition::applyBC() {
	int dim = face%3;
	int tag = 0;
	for (int i = 0; i < (int)ghostcells.size(); ++i) {
		for (int j = 0; j < (int)ghostcells[i].size(); ++j) {
			GridCell* ghost = ghostcells[i][j];
			GridCell* cptr = NULL;
			if (face < 3)
				cptr = ghostcells[i][j]->rjoin[dim]->rcell;
			else
				cptr = ghostcells[i][j]->ljoin[dim]->lcell;

			for (int k = 0; k < nghosts && cptr != NULL; ++k) {
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
	for (int i = 0; i < (int)ghostcells.size(); ++i) {
		for (int j = 0; j < (int)ghostcells[i].size(); ++j) {
			GridCell* ghost = ghostcells[i][j];
			while (ghost != NULL) {
				for(int iu = 0; iu < NU; ++iu) {
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
