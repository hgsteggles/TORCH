/* gridcell.cpp */

#include "gridcell.hpp"

#include <stddef.h>
#include <iostream>

GridJoin::GridJoin() {
	lcell = NULL;
	rcell = NULL;
	next = NULL;
	for(int i = 0; i < NU; ++i)
		F[i] = 0;
	for (int i = 0; i < 3; ++i)
		xj[i] = 0;
	area = 0;
	s_total++;
}
GridJoin::~GridJoin() {
	s_total--;
}
int GridJoin::s_total = 0;

GridCell::GridCell() : next(NULL), nextcausal(NULL), vol(0), ds(0), shellVol(0) {
	for(int dim = 0; dim < 3; ++dim){
		left[dim] = NULL;
		right[dim] = NULL;
		ljoin[dim] = NULL;
		rjoin[dim] = NULL;
	}
	for(int id = 0; id < 3; ++id){
		for(int iu = 0; iu < NU; ++iu){
			QL[id][iu] = 0;
			QR[id][iu] = 0;
		}
	}
	for(int i = 0; i < NU; ++i){
		U[i] = 0;
		Q[i] = 0;
		W[i] = 0;
	}
	for(int i = 0; i < NR; ++i){
		R[i] = 0;
	}
	for(int dim = 0; dim < 3; ++dim)
		xc[dim] = 0;
	for (int i = 0; i < 4; ++i) {
		NN[i] = NULL;
		NN_weights[i] = 0;
	}
	s_total++;
}
GridCell::~GridCell() {
	s_total--;
}
void GridCell::printInfo() {
	std::cout << "xc[0] = " << xc[0] << '\n';
	std::cout << "xc[1] = " << xc[1] << '\n';
	std::cout << "xc[2] = " << xc[2] << '\n';
	std::cout << "Q[iden] = " << Q[iden] << '\n';
	std::cout << "Q[ipre] = " << Q[ipre] << '\n';
	std::cout << "Q[ivel+0] = " << Q[ivel+0] << '\n';
	std::cout << "Q[ivel+1] = " << Q[ivel+1] << '\n';
	std::cout << "Q[ivel+2] = " << Q[ivel+2] << '\n';
	std::cout << "Q[ihii] = " << Q[ihii] << '\n';
	std::cout << "R[itau] = " << R[itau] << '\n';
	std::cout << "R[itauta] = " << R[itauta] << '\n';
	std::cout << "R[idtau] = " << R[idtau] << '\n';
	std::cout << "R[idtauta] = " << R[idtauta] << '\n';
	std::cout << "NN[0] = " << NN[0] << '\n';
	std::cout << "NN[1] = " << NN[1] << '\n';
	std::cout << "NN[2] = " << NN[2] << '\n';
	std::cout << "NN[3] = " << NN[3] << '\n';
	std::cout << "NN_weights[0] = " << NN_weights[0] << '\n';
	std::cout << "NN_weights[1] = " << NN_weights[1] << '\n';
	std::cout << "NN_weights[2] = " << NN_weights[2] << '\n';
	std::cout << "NN_weights[3] = " << NN_weights[3] << '\n';
	std::cout << "ljoin[0] = " << ljoin[0] << '\n';
	std::cout << "ljoin[1] = " << ljoin[1] << '\n';
	std::cout << "ljoin[2] = " << ljoin[2] << '\n';
	std::cout << "rjoin[0] = " << rjoin[0] << '\n';
	std::cout << "rjoin[1] = " << rjoin[1] << '\n';
	std::cout << "rjoin[2] = " << rjoin[2] << '\n';
	std::cout << "left[0] = " << left[0] << '\n';
	std::cout << "left[1] = " << left[1] << '\n';
	std::cout << "left[2] = " << left[2] << '\n';
	std::cout << "right[0] = " << right[0] << '\n';
	std::cout << "right[1] = " << right[1] << '\n';
	std::cout << "right[2] = " << right[2] << '\n';
}
int GridCell::s_total = 0;

/* SETS */
void GridCell::set_U(const int& index, const double& value) {U[index] = value;}
void GridCell::set_xcs(const int& x, const int& y, const int& z) {
	xc[0] = x;
	xc[1] = y;
	xc[2] = z;
}
/* GETS */
int GridCell::get_xc(const int& index) {return xc[index];}
double GridCell::get_U(const int& index) {return U[index];}
/* FUNCTIONS */
double GridCell::temperature() {
	double pre, rho;
	//molar_m = 1.0/(2.0*HIIfrac + (1.0-HIIfrac));
	pre = Q[ipre];
	rho = Q[iden];
	return pre/(GAS_CONST*rho);
}

