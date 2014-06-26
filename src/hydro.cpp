/**
 * @file hydro.cpp
 */

#include "hydro.h"
#include "grid3d.h"
#include "gridcell.h"
#include "boundary.h"
#include "radiation.h"
#include "parameters.h"
#include "slopelimiter.h"

#include "../Eigen/Dense"

#include <stddef.h> // NULL
#include <string.h>
#include <limits> // numeric_limits
#include <cmath> // sqrt, etc
#include <iostream> // cout, cerr

HydroDynamics::HydroDynamics(HydroParameters& hp, Grid3D& g3d) : hparams(hp), grid(g3d)
{
}

void HydroDynamics::globalWfromU() const {
	for(GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next){
		memcpy( (void*)cptr->W, (void*)cptr->U, NU * sizeof(double) );
		//for(int iu = 0; iu < NU; ++iu)
		//cptr->W[iu] = cptr->U[iu];
	}
}
void HydroDynamics::globalUfromW() const {
	for(GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next){
		memcpy( (void*)cptr->U, (void*)cptr->W, NU * sizeof(double) );
		//for(int iu = 0; iu < NU; ++iu)
		//cptr->U[iu] = cptr->W[iu];
	}
}
void HydroDynamics::globalQfromU() const {
	for(GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next)
		QfromU(cptr->Q, cptr->U);
}
void HydroDynamics::globalUfromQ() const {
	for(GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next)
		UfromQ(cptr->U, cptr->Q);
}
void HydroDynamics::reconstruct() const {
	double dQdr1[NU], dQdr2[NU], dQdr[NU];
	for(GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next){
		for(int dim = 0; dim < grid.gparams->ND; ++dim){
			for(int iq = 0; iq < NU; ++iq){
				if(cptr->left[dim] != NULL)
					dQdr1[iq] = (cptr->Q[iq] - cptr->left[dim]->Q[iq])/grid.dx[dim];
				else
					dQdr1[iq] = (cptr->Q[iq] - cptr->ljoin[dim]->lcell->Q[iq])/grid.dx[dim];
				if(cptr->right[dim] != NULL)
					dQdr2[iq] = (cptr->right[dim]->Q[iq] - cptr->Q[iq])/grid.dx[dim];
				else
					dQdr2[iq] = (cptr->rjoin[dim]->rcell->Q[iq] - cptr->Q[iq])/grid.dx[dim];
				dQdr[iq] = SlopeLimiter::SF(dQdr1[iq], dQdr2[iq]);
				cptr->QL[dim][iq] = cptr->Q[iq] - 0.5*(grid.dx[dim])*dQdr[iq];
				cptr->QR[dim][iq] = cptr->Q[iq] + 0.5*(grid.dx[dim])*dQdr[iq];
			}
			//UfromQ(cptr->UL[dim], cptr->QL[dim]);
			//UfromQ(cptr->UR[dim], cptr->QR[dim]);
		}
	}
	for (int ib = 0; ib < (int)grid.boundaries.size(); ++ib) {
		int face = grid.boundaries[ib]->face;
		int dim = face%3;
		for (int i = 0; i < (int)grid.boundaries[ib]->ghostcells.size(); ++i) {
			for (int j = 0; j < (int)grid.boundaries[ib]->ghostcells[i].size(); ++j) {
				GridCell* bcptr = grid.boundaries[ib]->ghostcells[i][j];
				GridCell* cptr = NULL;
				if (face < 3) {
					cptr = bcptr->rjoin[dim]->rcell;
					for(int iq = 0; iq < NU; ++iq){
						dQdr1[iq] = (bcptr->Q[iq] - bcptr->left[dim]->Q[iq])/grid.dx[dim];
						dQdr2[iq] = (cptr->Q[iq] - bcptr->Q[iq])/grid.dx[dim];
					}
				}
				else {
					cptr = bcptr->ljoin[dim]->lcell;
					for (int iq = 0; iq < NU; ++iq) {
						dQdr1[iq] = (bcptr->Q[iq] - cptr->Q[iq])/grid.dx[dim];
						dQdr2[iq] = (bcptr->right[dim]->Q[iq] - bcptr->Q[iq])/grid.dx[dim];
					}
				}
				for (int iq = 0; iq < NU; ++iq) {
					dQdr[iq] = SlopeLimiter::SF(dQdr1[iq], dQdr2[iq]);
					bcptr->QL[dim][iq] = bcptr->Q[iq] - 0.5*(grid.dx[dim])*dQdr[iq];
					bcptr->QR[dim][iq] = bcptr->Q[iq] + 0.5*(grid.dx[dim])*dQdr[iq];
				}

				//UfromQ(bcptr->UL[dim], bcptr->QL[dim]);
				//UfromQ(bcptr->UR[dim], bcptr->QR[dim]);
			}
		}
	}
}
double HydroDynamics::av(const double& a, const double& b) const {
	double result = 0;
	double nom = a*a*b + a*b*b;
	double denom = a*a + b*b;
	if(denom > std::numeric_limits<double>::min())
		result = nom/denom;
	return result;
}

void HydroDynamics::UfromQ(double u[], double q[]) const {
	double ke = 0;
	u[iden] = q[iden] = std::max(q[iden], hparams.DFLOOR);
	for(int dim = 0; dim < grid.gparams->ND; ++dim){
		u[ivel+dim] = q[ivel+dim]*q[iden];
		ke += 0.5*q[iden]*q[ivel+dim]*q[ivel+dim];
	}
	q[ipre] = std::max(q[ipre], hparams.PFLOOR);
	u[ipre] = q[ipre]/(hparams.GAMMA - 1.0) + ke;
	q[ihii] = std::max(std::min(q[ihii], 1.0), 0.0);
	u[ihii] = q[ihii]*q[iden];
}
void HydroDynamics::QfromU(double q[], double u[]) const {
	q[iden] = u[iden] = std::max(u[iden], hparams.DFLOOR);
	double ke = 0;
	for (int dim = 0; dim < grid.gparams->ND; ++dim) {
		q[ivel+dim] = u[ivel+dim]/u[iden];
		ke += 0.5*u[ivel+dim]*u[ivel+dim]/u[iden];
	}
	u[ipre] = std::max(u[ipre], hparams.PFLOOR/(hparams.GAMMA - 1.0) + ke);
	q[ipre] = (u[ipre] - ke)*(hparams.GAMMA - 1.0);
	u[ihii] = std::max(std::min(u[ihii], u[iden]), 0.0);
	q[ihii] = u[ihii]/u[iden];
}
void HydroDynamics::FfromU(double f[], double u[], const int& dim) const {
	double ke = 0, pressure;
	for(int id = 0; id < grid.gparams->ND; ++id)
		ke += 0.5*u[ivel+id]*u[ivel+id]/u[iden];
	pressure = (u[ipre] - ke)*(hparams.GAMMA - 1.0);
	f[iden] = u[ivel+dim];
	f[ivel+dim] = u[ivel+dim]*u[ivel+dim]/u[iden] + pressure;
	f[ivel+(dim+1)%3] = u[ivel+(dim+1)%3]*u[ivel+dim]/u[iden];
	f[ivel+(dim+2)%3] = u[ivel+(dim+2)%3]*u[ivel+dim]/u[iden];
	f[ipre] = u[ivel+dim]*(u[ipre] + pressure)/u[iden];
	f[ihii] = u[ivel+dim]*u[ihii]/u[iden];
}
double HydroDynamics::soundSpeed(const double& pre, const double& den) const {
	return std::sqrt(hparams.GAMMA*pre/den);
}
double HydroDynamics::CFL(const double& dt_max) const {
	double tmin = dt_max;
	for(GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next){
		double inv_t = 0;
		double ss = soundSpeed(cptr->Q[ipre], cptr->Q[iden]);
		for(int dim = 0; dim < grid.gparams->ND; ++dim)
			inv_t += (fabs(cptr->Q[ivel+dim]) + ss)/grid.dx[dim];
		if(inv_t == 0)
			inv_t = 1.0/dt_max;
		tmin = std::min(tmin, 0.5/inv_t);
		if (0.5/inv_t < 1.0e-8) {
			std::cout << "ERROR: hydro step too small. Pressure = " << cptr->Q[ipre];
			std::cout << ", Density = " << cptr->Q[iden] << ", Cs = " << ss;
			std::cout << ", v = {" << cptr->Q[ivel+0] << ", " << cptr->Q[ivel+1] << ", " << cptr->Q[ivel+2] << "}\n";
		}
		/*
		if (1.0/inv_t > dt_max) {
			std::cout << "ERROR: hydro step too large. Pressure = " << cptr->Q[ipre];
			std::cout << ", Density = " << cptr->Q[iden] << ", Cs = " << ss;
			std::cout << ", v = {" << cptr->Q[ivel+0] << ", " << cptr->Q[ivel+1] << ", " << cptr->Q[ivel+2] << "}";
			std::cout << ", dx = {" << grid.dx[0] << ", " << grid.dx[1] << ", " << grid.dx[2] << "}";
			std::cout << ", gamma = " << hparams.GAMMA << "\n";
		}
		*/
	}
	return std::min(tmin, dt_max);
}

void HydroDynamics::calcFluxes(const int& order) const {
	if(order == 1)
		reconstruct();
	for(int dim = 0; dim < grid.gparams->ND; ++dim){
		for(GridJoin* jptr = grid.fjoin[dim]; jptr != NULL; jptr = jptr->next){
			if(order == 1)
				RotatedHLLC(jptr->lcell->QR[dim], jptr->F, jptr->rcell->QL[dim], dim);
			else if (order == 0)
				RotatedHLLC(jptr->lcell->Q, jptr->F, jptr->rcell->Q, dim);
			else {
				std::cerr << "ERROR: invalid order(=" << order << ") passed into calcFluxes. Valid orders = {0, 1}." << std::endl;
				exit(EXIT_FAILURE);
			}
		}
	}
}

void HydroDynamics::advSolution(const double& dt) {
	for (GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next) {
		//Fluxes.
		for (int dim = 0; dim < grid.gparams->ND; ++dim) {
			for (int i = 0; i < NU; ++i) {
				if(cptr->ljoin[dim] != NULL && cptr->rjoin[dim] != NULL) {
					cptr->U[i] += (dt*cptr->ljoin[dim]->area/cptr->vol)*cptr->ljoin[dim]->F[i];
					cptr->U[i] -= (dt*cptr->rjoin[dim]->area/cptr->vol)*cptr->rjoin[dim]->F[i];
				}
				else {
					std::cerr << "ERROR: ljoin[dim] or rjoin[dim] is NULL in Hydrodynamics::advSolution()." << '\n';
					exit(EXIT_FAILURE);
				}
			}
		}
	}
}

void HydroDynamics::updateSrcTerms() const {
	bool rtp = (grid.gparams->GEOMETRY == SPHERICAL);
	bool rzp = (grid.gparams->GEOMETRY== CYLINDRICAL);

	for (GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next) {
		//Fluxes.
		for (int dim = 0; dim < grid.gparams->ND; ++dim) {
			for (int i = 0; i < NU; ++i) {
				if(cptr->ljoin[dim] != NULL && cptr->rjoin[dim] != NULL) {
					cptr->UDOT[i] += (cptr->ljoin[dim]->area/cptr->vol)*cptr->ljoin[dim]->F[i];
					cptr->UDOT[i] -= (cptr->rjoin[dim]->area/cptr->vol)*cptr->rjoin[dim]->F[i];
				}
				else {
					std::cerr << "ERROR: ljoin[dim] or rjoin[dim] is NULL in Hydrodynamics::advSolution()." << '\n';
					exit(EXIT_FAILURE);
				}
			}
		}
		//Geometric.
		if (rzp) {
			double r = grid.dx[0]*(cptr->xc[0] + 0.5);
			cptr->UDOT[ivel+0] += cptr->Q[ipre]/r;
		}
		else if (rtp){
			double r, area1, area2;
			area1 = cptr->rjoin[0]->area;
			area2 = cptr->ljoin[0]->area;
			r = cptr->vol/(area1 - area2);
			cptr->UDOT[ivel+0] += cptr->Q[ipre]/r;
		}
	}
}

void HydroDynamics::updateBoundaries() const {
	for(int i = 0; i < (int)grid.boundaries.size(); ++i)
		grid.boundaries[i]->applyBC();
}

void HydroDynamics::HLLC(double Q_l[], double F[], double Q_r[], const int& dim) const {
	double U_cl[NU], U_cr[NU], U_l[NU], U_r[NU], F_l[NU], F_r[NU], F_cl[NU], F_cr[NU];
	for(int i = 0; i < NU; ++i)
		U_cl[i] = U_cr[i] = U_l[i] = U_r[i] = F_l[i] = F_r[i] = F_cl[i] = F_cr[i] = 0;
	UfromQ(U_l, Q_l);
	UfromQ(U_r, Q_r);
	FfromU(F_l, U_l, dim);
	FfromU(F_r, U_r, dim);

	/**
	 * Davis Estimates
	 * Not recommended for practical computations
	 * Toro: Riemann Solvers + Numerical Methods for Fluid Dynamics (pg 328)
	 * double S_l = Q_l[ivel+dim] - soundSpeed(Q_l[ipre], Q_l[iden]);
	 * double S_r = Q_r[ivel+dim] + soundSpeed(Q_r[ipre], Q_r[iden]);
	 */

	/**
	 * Einfeldt Estimates
	 * Used in HLLE solver
	 * Toro: Riemann Solvers + Numerical Methods for Fluid Dynamics (pg 328)
	 */
	double sqrtrho_l = std::sqrt(Q_l[iden]);
	double sqrtrho_r = std::sqrt(Q_r[iden]);
	double u_tilde = (sqrtrho_l*Q_l[ivel+dim] + sqrtrho_r*Q_r[ivel+dim])/(sqrtrho_l + sqrtrho_r);
	double a_l = soundSpeed(Q_l[ipre], Q_l[iden]);
	double a_r = soundSpeed(Q_r[ipre], Q_r[iden]);
	double nu2 = 0.5*sqrtrho_r*sqrtrho_l/((sqrtrho_l + sqrtrho_r)*(sqrtrho_l + sqrtrho_r));
	double dsqrd = (sqrtrho_l*a_l*a_l + sqrtrho_r*a_r*a_r)/(sqrtrho_l + sqrtrho_r);
	dsqrd += nu2*(Q_r[ivel+dim] - Q_l[ivel+dim])*(Q_r[ivel+dim] - Q_l[ivel+dim]);
	double d = std::sqrt(dsqrd);
	double S_l = u_tilde - d;
	double S_r = u_tilde + d;
	/****************************/

	// S_c = (P_r-P_l+D_l*u_l*(S_l-u_l) - D_r*u_r*(S_r-u_r))/(D_l*(S_l-u_l) - D_r*(S_r-u_r));
	double S_c = (Q_r[ipre]-Q_l[ipre]+Q_l[iden]*Q_l[ivel+dim]*(S_l-Q_l[ivel+dim])-Q_r[iden]*Q_r[ivel+dim]*(S_r-Q_r[ivel+dim]))
					/ (Q_l[iden]*(S_l-Q_l[ivel+dim])-Q_r[iden]*(S_r-Q_r[ivel+dim]));
	double A_r = Q_r[iden]*(S_r-Q_r[ivel+dim])/(S_r-S_c);
	U_cr[iden] = A_r;
	for(int id = 0; id < grid.gparams->ND; ++id)
		U_cr[ivel+id] = A_r*Q_r[ivel+id];
	U_cr[ivel+dim] = A_r*S_c;
	U_cr[ipre] = A_r*((U_r[ipre]/Q_r[iden]) + (S_c-Q_r[ivel+dim])*(S_c + Q_r[ipre]/(Q_r[iden]*(S_r-Q_r[ivel+dim]))));
	U_cr[ihii] = A_r*Q_r[ihii];

	double A_l = Q_l[iden]*(S_l-Q_l[ivel+dim])/(S_l-S_c);
	U_cl[iden] = A_l;
	for(int id = 0; id < grid.gparams->ND; ++id)
		U_cl[ivel+id] = A_l*Q_l[ivel+id];
	U_cl[ivel+dim] = A_l*S_c;
	U_cl[ipre] = A_l*((U_l[ipre]/Q_l[iden]) + (S_c-Q_l[ivel+dim])*(S_c+Q_l[ipre]/(Q_l[iden]*(S_l-Q_l[ivel+dim]))));
	U_cl[ihii] = A_l*Q_l[ihii];
	for(int i = 0; i < NU; ++i){
		F_cl[i] = F_l[i]+S_l*(U_cl[i]-U_l[i]);
		F_cr[i] = F_r[i]+S_r*(U_cr[i]-U_r[i]);
	}
	if(S_l >= 0){
		for(int i = 0; i < NU; ++i)
			F[i] = F_l[i];
	}
	else if(S_l <= 0 && S_c >= 0){
		for(int i = 0; i < NU; ++i)
			F[i] = F_cl[i];
	}
	else if (S_c <= 0 && S_r >= 0){
		for(int i = 0; i < NU; ++i)
			F[i] = F_cr[i];
	}
	else if(S_r <= 0){
		for(int i = 0; i < NU; ++i)
			F[i] = F_r[i];
	}
	if(std::isnan(F[ihii])) {
		std::cerr << '\n';
		for(int i = 0; i < NU; ++i){
			std::cerr << "Q_l[" << i << "] = " << Q_l[i] << '\n';
			std::cerr << "Q_r[" << i << "] = " << Q_r[i] << '\n';
		}
		for(int i = 0; i < NU; ++i){
			std::cerr << "U_l[" << i << "] = " << U_l[i] << '\n';
			std::cerr << "U_r[" << i << "] = " << U_r[i] << '\n';
		}
		for(int i = 0; i < NU; ++i){
			std::cerr << "F_l[" << i << "] = " << F_l[i] << '\n';
			std::cerr << "F_r[" << i << "] = " << F_r[i] << '\n';
		}
		for(int i = 0; i < NU; ++i){
			std::cerr << "U_cl[" << i << "] = " << U_cl[i] << '\n';
			std::cerr << "U_cr[" << i << "] = " << U_cr[i] << '\n';
		}
		std::cerr << "a_l = " << soundSpeed(Q_l[ipre], Q_l[iden]) << '\n';
		std::cerr << "a_r = " << soundSpeed(Q_r[ipre], Q_r[iden]) << '\n';
		std::cerr << "S_l = " << S_l << '\n';
		std::cerr << "S_r = " << S_r << '\n';
		std::cerr << "S_c = " << S_c << '\n';
		for(int iu = 0; iu < NU; ++iu)
			std::cerr << "F[" << iu << "] = " << F[iu] << '\n';
	}
}

void HydroDynamics::HLL(double Q_l[], double F[], double Q_r[], const int& dim) const {
	double U_l[NU], U_r[NU], F_l[NU], F_r[NU], F_c[NU];
	for(int i = 0; i < NU; ++i)
		U_l[i] = U_r[i] = F_l[i] = F_r[i] = F_c[i] = 0;
	UfromQ(U_l, Q_l);
	UfromQ(U_r, Q_r);
	FfromU(F_l, U_l, dim);
	FfromU(F_r, U_r, dim);
	/**
	 * Einfeldt Estimates
	 * Used in HLLE solver
	 * Toro: Riemann Solvers + Numerical Methods for Fluid Dynamics (pg 328)
	 */
	double sqrtrho_l = std::sqrt(Q_l[iden]);
	double sqrtrho_r = std::sqrt(Q_r[iden]);
	double u_tilde = (sqrtrho_l*Q_l[ivel+dim] + sqrtrho_r*Q_r[ivel+dim])/(sqrtrho_l + sqrtrho_r);
	double a_l = soundSpeed(Q_l[ipre], Q_l[iden]);
	double a_r = soundSpeed(Q_r[ipre], Q_r[iden]);
	double nu2 = 0.5*sqrtrho_r*sqrtrho_l/((sqrtrho_l + sqrtrho_r)*(sqrtrho_l + sqrtrho_r));
	double dsqrd = (sqrtrho_l*a_l*a_l + sqrtrho_r*a_r*a_r)/(sqrtrho_l + sqrtrho_r);
	dsqrd += nu2*(Q_r[ivel+dim] - Q_l[ivel+dim])*(Q_r[ivel+dim] - Q_l[ivel+dim]);
	double d = std::sqrt(dsqrd);
	double S_l = u_tilde - d;
	double S_r = u_tilde + d;

	for(int i = 0; i < NU; ++i)
		F_c[i] = (S_r*F_l[i] - S_l*F_r[i] + S_l*S_r*(U_r[i]-U_l[i]))/(S_r-S_l);
	if(S_l >= 0){
		for(int i = 0; i < NU; ++i)
			F[i] = F_l[i];
	}
	else if(S_l <= 0 && S_r >= 0){
		for(int i = 0; i < NU; ++i)
			F[i] = F_c[i];
	}
	else if(S_r <= 0){
		for(int i = 0; i < NU; ++i)
			F[i] = F_r[i];
	}
}

Eigen::Matrix<double, 3, 3> getRotationMatrix(Eigen::Matrix<double, 3, 1> l, double cos) {
	Eigen::Matrix<double, 3, 3> R;
	double sin = std::sqrt(1.0 - cos*cos);
	R << l(0)*l(0) + (l(1)*l(1) + l(2)*l(2))*cos, l(0)*l(1)*(1.0-cos) - l(2)*sin         , l(0)*l(2)*(1.0-cos) + l(1)*sin         ,
			l(0)*l(1)*(1.0-cos) + l(2)*sin         , l(1)*l(1) + (l(0)*l(0) + l(2)*l(2))*cos, l(1)*l(2)*(1.0-cos) - l(0)*sin         ,
			l(0)*l(2)*(1.0-cos) - l(1)*sin         , l(1)*l(2)*(1.0-cos) + l(0)*sin         , l(2)*l(2) + (l(0)*l(0) + l(1)*l(1))*cos;
	return R;
}
Eigen::Matrix<double, 3, 3> getInverseRotationMatrix(Eigen::Matrix<double, 3, 1> l, double cos) {
	Eigen::Matrix<double, 3, 3> R;
	double sin = -std::sqrt(1.0 - cos*cos);
	R << l(0)*l(0) + (l(1)*l(1) + l(2)*l(2))*cos, l(0)*l(1)*(1.0-cos) - l(2)*sin         , l(0)*l(2)*(1.0-cos) + l(1)*sin         ,
			l(0)*l(1)*(1.0-cos) + l(2)*sin         , l(1)*l(1) + (l(0)*l(0) + l(2)*l(2))*cos, l(1)*l(2)*(1.0-cos) - l(0)*sin         ,
			l(0)*l(2)*(1.0-cos) - l(1)*sin         , l(1)*l(2)*(1.0-cos) + l(0)*sin         , l(2)*l(2) + (l(0)*l(0) + l(1)*l(1))*cos;
	return R;
}

void rotate(Eigen::Matrix<double, 3, 3> R, double v[]) {
	Eigen::Matrix<double, 3, 1> vel(v[ivel+0], v[ivel+1], v[ivel+2]);
	vel = R*vel;
	for (int i = 0; i < 3; ++i)
		v[ivel+i] = vel(i);
}

void HydroDynamics::RotatedHLLC(double Q_l[], double F[], double Q_r[], const int& dim) const {
	bool debugger = false;
	bool pureHLLC = false, pureHLL = false;
	Eigen::Matrix<double, 3, 1> d(dim==0 ? 1 : 0, dim==1 ? 1 : 0, dim==2 ? 1 : 0);
	Eigen::Matrix<double, 3, 1> n1(Q_r[ivel+0]-Q_l[ivel+0], Q_r[ivel+1]-Q_l[ivel+1], Q_r[ivel+2]-Q_l[ivel+2]);
	Eigen::Matrix<double, 3, 1> n2, axis1, axis2;
	if (!(Q_r[ivel+0] == Q_r[ivel+0])) {
		std::cerr << "ERROR: Q_r is nan in RotatedHLLC." << std::endl;
		exit(EXIT_FAILURE);
	}
	if (n1(dim) < 0)
		n1 = -1*n1;
	if (n1.norm() < 1.0e-6)
		pureHLLC = true;
	else {
		n1.normalize();
		axis1 = (n1.cross(d));
		if (axis1.norm() == 0)
			pureHLL = true;
		else {
			n2 = (n1.cross(d)).cross(n1);
			if (n2.norm() == 0)
				pureHLL = true;
			else {
				n2.normalize();
				axis2 = (n2.cross(d));
				if (axis2.norm() == 0)
					pureHLLC = true;
				else
					axis2.normalize();
			}
		}
	}
	if (!pureHLLC && !pureHLL) {
		double F1[NU], F2[NU], alpha1=std::abs(d.dot(n1)), alpha2=std::abs(d.dot(n2));
		if (axis1.norm() != 0)
			axis1.normalize();
		else
			std::cerr << "ERROR: axis1 is a zero vector." << std::endl;
		Eigen::Matrix<double, 3, 3> R1 = getRotationMatrix(axis1, alpha1);
		rotate(R1, Q_l);
		rotate(R1, Q_r);
		HLL(Q_l, F1, Q_r, dim);
		Eigen::Matrix<double, 3, 3> R1_inv = getInverseRotationMatrix(axis1, alpha1);
		rotate(R1_inv, Q_l);
		rotate(R1_inv, Q_r);
		rotate(R1_inv, F1);
		Eigen::Matrix<double, 3, 3> R2 = getRotationMatrix(axis2, alpha2);
		rotate(R2, Q_l);
		rotate(R2, Q_r);
		HLLC(Q_l, F2, Q_r, dim);
		Eigen::Matrix<double, 3, 3> R2_inv = getInverseRotationMatrix(axis2, alpha2);
		rotate(R2_inv, Q_l);
		rotate(R2_inv, Q_r);
		rotate(R2_inv, F2);
		for (int iu = 0; iu < NU; ++iu)
			F[iu] = alpha1*F1[iu] + alpha2*F2[iu];
		if (debugger) {
			std::cerr << "F1 = [ " << F1[0] << ", " << F1[1] << ", " << F1[2] << ", " << F1[3] << ", " << F1[4] << ", " << F1[5] << " ]\n";
			std::cerr << "F2 = [ " << F2[0] << ", " << F2[1] << ", " << F2[2] << ", " << F2[3] << ", " << F2[4] << ", " << F2[5] << " ]\n";
		}
		if (std::isnan(F1[1]) && debugger) {
			std::cerr << "d = " << d << std::endl;
			std::cerr << "n1 = " << n1 << std::endl;
			std::cerr << "n2 = " << n2 << std::endl;
		}
	}
	else if (pureHLL) {
		HLL(Q_l, F, Q_r, dim);
		if(std::isnan(F[ihii]) && debugger) {
			std::cerr << "ERROR: nan detected after call to HLL. Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	else if (pureHLLC) {
		HLLC(Q_l, F, Q_r, dim);
		if(std::isnan(F[ihii]) && debugger) {
			std::cerr << "ERROR: nan detected after call to HLLC. Exiting..." << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

void HydroDynamics::applySrcTerms(const double& dt) const {
	/* Geometry */
	bool rtp = (grid.gparams->GEOMETRY == SPHERICAL);
	bool rzp = (grid.gparams->GEOMETRY== CYLINDRICAL);
	//bool xyz = (grid.gparams->GEOMETRY == CARTESIAN);
	for (GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next) {
		/* Geometric source term (if any) */
		if (rzp) {
			double r = grid.dx[0]*(cptr->xc[0] + 0.5);
			cptr->U[ivel+0] += dt*cptr->Q[ipre]/r;
		}
		else if (rtp){
			double r, area1, area2;
			area1 = cptr->rjoin[0]->area;
			area2 = cptr->ljoin[0]->area;
			r = cptr->vol/(area1 - area2);
			cptr->U[ivel+0] += dt*cptr->Q[ipre]/r;
		}
	}
}
void HydroDynamics::fixSolution() const {
	for (GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next) {
		cptr->U[iden] = std::max(cptr->U[iden], hparams.DFLOOR);
		double ke = 0.0;
		for(int dim = 0; dim < grid.gparams->ND; ++dim)
			ke += 0.5*cptr->U[ivel+dim]*cptr->U[ivel+dim]/cptr->U[iden];
		cptr->U[ipre] = std::max(cptr->U[ipre], hparams.PFLOOR/(hparams.GAMMA - 1.0) + ke);
		cptr->U[ihii] = std::max(std::min(cptr->U[ihii], cptr->U[iden]), 0.0);
	}
}
void HydroDynamics::Qisnan(const int& id, const int& i, const int& xc, const int& yc, const int& zc) const {
	GridCell* cptr = grid.locate(xc, yc, zc);
	if(std::isnan(cptr->Q[i]))
		std::cerr << id << '\n';
}
