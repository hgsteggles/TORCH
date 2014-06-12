/**
 * @file gridcell.cpp
 */

#include "gridcell.hpp"

#include <stddef.h>
#include <iostream>

/**
 * @brief The default GridJoin constructor.
 * Provides all attributes with safe values.
 */
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

/**
 * @brief The Gridjoin destructor.
 */
GridJoin::~GridJoin() {
	s_total--;
}

int GridJoin::s_total = 0;

/**
 * @brief Default GridCell constructor.
 * Provides all attributes with safe values.
 */
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
		UDOT[i] = 0;
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

/**
 * @brief GridCell destructor.
 * Does NOT delete objects that this GridCell object points to.
 */
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

/**
 * @brief Setter for GridCell::U.
 * @param index
 * @param value
 */
void GridCell::set_U(const int& index, const double& value) {U[index] = value;}

/**
 * @brief Setter for GridCell:xc.
 * @param x The x grid coordinate.
 * @param y The y grid coordinate.
 * @param z The z grid coordinate.
 */
void GridCell::set_xcs(const int& x, const int& y, const int& z) {
	xc[0] = x;
	xc[1] = y;
	xc[2] = z;
}

/**
 * @brief Getter for GridCell::xc.
 * @param i The grid coordinate to be returned.
 * @return The location of this GridCell object on grid coordinate i.
 */
int GridCell::get_xc(const int& index) {return xc[index];}

/**
 * @brief Getter for GridCell::U.
 * @param index The index for the fluid variable to be returned.
 * @return The value of the fluid variable.
 */
double GridCell::get_U(const int& index) {return U[index];}

/**
 * @brief Returns the temperature of this GridCell object.
 * @return Temperature.
 */
double GridCell::temperature() {
	double pre, rho;
	//molar_m = 1.0/(2.0*HIIfrac + (1.0-HIIfrac));
	pre = Q[ipre];
	rho = Q[iden];
	return pre/(GAS_CONST*rho);
}

