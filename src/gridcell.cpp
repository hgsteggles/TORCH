#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <math.h>
#include "gridcell.hpp"
#include <map>

using namespace std;

GridJoin::GridJoin() {
	lcell = NULL;
	rcell = NULL;
	next = NULL;
	for(int i = 0; i < NU; i++)
		F[i] = 0;
	for (int i = 0; i < 3; i++)
		xj[i] = 0;
	area = 0;
	s_total++;
}
GridJoin::~GridJoin() {
	s_total--;
}
int GridJoin::s_total = 0;

GridCell::GridCell() {
	for(int dim = 0; dim < 3; dim++){
		left[dim] = NULL;
		right[dim] = NULL;
		ljoin[dim] = NULL;
		rjoin[dim] = NULL;
	}
	next = NULL;
	for(int id = 0; id < 3; id++){
		for(int iu = 0; iu < NU; iu++){
			UL[id][iu] = 0;
			UR[id][iu] = 0;
			QL[id][iu] = 0;
			QR[id][iu] = 0;
		}
	}
	for(int i = 0; i < NU; i++){
		U[i] = 0;
		Q[i] = 0;
		W[i] = 0;
	}
	for(int i = 0; i < NR; i++){
		R[i] = 0;
	}
	for(int dim = 0; dim < 3; dim++)
		xc[dim] = 0;
	vol = 0;
	s_total++;
}
GridCell::~GridCell() {
	s_total--;
}
int GridCell::s_total = 0;

/* SETS */
void GridCell::set_U(int index, double value) {U[index] = value;}
void GridCell::set_xcs(int x, int y, int z) {
	xc[0] = x;
	xc[1] = y;
	xc[2] = z;
}
/* GETS */
int GridCell::get_xc(int index) {return xc[index];}
double GridCell::get_U(int index) {return U[index];}
/* FUNCTIONS */
double GridCell::temperature() {
	double pre, rho;
	//molar_m = 1.0/(2.0*HIIfrac + (1.0-HIIfrac));
	pre = Q[ipre];
	rho = Q[iden];
	return pre/(GAS_CONST*rho);
}

