/* hydro.cpp */

#include "hydro.hpp"
#include "grid3d.hpp"
#include "gridcell.hpp"
#include "boundary.hpp"
#include "rtmodule.hpp"
#include "parameters.hpp"

#include <stddef.h> // NULL
#include <string.h>
#include <limits> // numeric_limits
#include <cmath> // sqrt, etc
#include <iostream> // cout, cerr

HydroDynamics::HydroDynamics(const HydroParameters& hp, Grid3D* grid) : gptr(grid) {
	GAMMA = hp.GAMMA;
	DFLOOR = hp.DFLOOR;
	PFLOOR = hp.PFLOOR;
	DTMAX = hp.DTMAX;
}

void HydroDynamics::globalWfromU() const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
		memcpy( (void*)cptr->W, (void*)cptr->U, NU * sizeof(double) );
		//for(int iu = 0; iu < NU; ++iu)
			//cptr->W[iu] = cptr->U[iu];
	}
}
void HydroDynamics::globalUfromW() const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
		memcpy( (void*)cptr->U, (void*)cptr->W, NU * sizeof(double) );
		//for(int iu = 0; iu < NU; ++iu)
			//cptr->U[iu] = cptr->W[iu];
	}
}
void HydroDynamics::globalQfromU() const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next)
		QfromU(cptr->Q, cptr->U);
}
void HydroDynamics::globalUfromQ() const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next)
		UfromQ(cptr->U, cptr->Q);
}
void HydroDynamics::reconstruct() const {
	double dQdr1[NU], dQdr2[NU], dQdr[NU];
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
		for(int dim = 0; dim < gptr->ND; ++dim){
			for(int iq = 0; iq < NU; ++iq){
				if(cptr->left[dim] != NULL)
					dQdr1[iq] = (cptr->Q[iq] - cptr->left[dim]->Q[iq])/gptr->dx[dim];
				else
					dQdr1[iq] = (cptr->Q[iq] - cptr->ljoin[dim]->lcell->Q[iq])/gptr->dx[dim];
				if(cptr->right[dim] != NULL)
					dQdr2[iq] = (cptr->right[dim]->Q[iq] - cptr->Q[iq])/gptr->dx[dim];
				else
					dQdr2[iq] = (cptr->rjoin[dim]->rcell->Q[iq] - cptr->Q[iq])/gptr->dx[dim];
				dQdr[iq] = av(dQdr1[iq], dQdr2[iq]);
				cptr->QL[dim][iq] = cptr->Q[iq] - 0.5*(gptr->dx[dim])*dQdr[iq];
				cptr->QR[dim][iq] = cptr->Q[iq] + 0.5*(gptr->dx[dim])*dQdr[iq];
			}
			//UfromQ(cptr->UL[dim], cptr->QL[dim]);
			//UfromQ(cptr->UR[dim], cptr->QR[dim]);
		}
	}
	for (int ib = 0; ib < (int)gptr->boundaries.size(); ++ib) {
		int face = gptr->boundaries[ib]->face;
		int dim = face%3;
		for (int i = 0; i < (int)gptr->boundaries[ib]->ghostcells.size(); ++i) {
			for (int j = 0; j < (int)gptr->boundaries[ib]->ghostcells[i].size(); ++j) {
				GridCell* bcptr = gptr->boundaries[ib]->ghostcells[i][j];
				GridCell* cptr = NULL;
				if (face < 3) {
					cptr = bcptr->rjoin[dim]->rcell;
					for(int iq = 0; iq < NU; ++iq){
						dQdr1[iq] = (bcptr->Q[iq] - bcptr->left[dim]->Q[iq])/gptr->dx[dim];
						dQdr2[iq] = (cptr->Q[iq] - bcptr->Q[iq])/gptr->dx[dim];
					}
				}
				else {
					cptr = bcptr->ljoin[dim]->lcell;
					for (int iq = 0; iq < NU; ++iq) {
						dQdr1[iq] = (bcptr->Q[iq] - cptr->Q[iq])/gptr->dx[dim];
						dQdr2[iq] = (bcptr->right[dim]->Q[iq] - bcptr->Q[iq])/gptr->dx[dim];
					}
				}
				for (int iq = 0; iq < NU; ++iq) {
					dQdr[iq] = av(dQdr1[iq], dQdr2[iq]);
					bcptr->QL[dim][iq] = bcptr->Q[iq] - 0.5*(gptr->dx[dim])*dQdr[iq];
					bcptr->QR[dim][iq] = bcptr->Q[iq] + 0.5*(gptr->dx[dim])*dQdr[iq];
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
	u[iden] = q[iden] = std::max(q[iden], DFLOOR);
	for(int dim = 0; dim < gptr->ND; ++dim){
		u[ivel+dim] = q[ivel+dim]*q[iden];
		ke += 0.5*q[iden]*q[ivel+dim]*q[ivel+dim];
	}
	q[ipre] = std::max(q[ipre], PFLOOR);
	u[ipre] = q[ipre]/(GAMMA - 1.0) + ke;
	q[ihii] = std::max(std::min(q[ihii], 1.0), 0.0);
	u[ihii] = q[ihii]*q[iden];
}
void HydroDynamics::QfromU(double q[], double u[]) const {
	q[iden] = u[iden] = std::max(u[iden], DFLOOR);
	double ke = 0;
	for(int dim = 0; dim < gptr->ND; ++dim){
		q[ivel+dim] = u[ivel+dim]/u[iden];
		ke += 0.5*u[ivel+dim]*u[ivel+dim]/u[iden];
	}
	u[ipre] = std::max(u[ipre], PFLOOR/(GAMMA - 1.0) + ke);
	q[ipre] = (u[ipre] - ke)*(GAMMA - 1.0);
	u[ihii] = std::max(std::min(u[ihii], u[iden]), 0.0);
	q[ihii] = u[ihii]/u[iden];
}
void HydroDynamics::FfromU(double f[], double u[], const int& dim) const {
	double ke = 0, pressure;
	for(int id = 0; id < gptr->ND; ++id)
		ke += 0.5*u[ivel+id]*u[ivel+id]/u[iden];
	pressure = (u[ipre] - ke)*(GAMMA - 1.0);
	f[iden] = u[ivel+dim];
	f[ivel+dim] = u[ivel+dim]*u[ivel+dim]/u[iden] + pressure;
	f[ivel+(dim+1)%3] = u[ivel+(dim+1)%3]*u[ivel+dim]/u[iden];
	f[ivel+(dim+2)%3] = u[ivel+(dim+2)%3]*u[ivel+dim]/u[iden];
	f[ipre] = u[ivel+dim]*(u[ipre] + pressure)/u[iden];
	f[ihii] = u[ivel+dim]*u[ihii]/u[iden];
}
double HydroDynamics::soundSpeed(const double& pre, const double& den) const {
	return sqrt(GAMMA*pre/den);
}
double HydroDynamics::CFL() const {
	double tmin = DTMAX;
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
		double inv_t = 0;
		double ss = soundSpeed(cptr->Q[ipre], cptr->Q[iden]);
		for(int dim = 0; dim < gptr->ND; ++dim)
			inv_t += (fabs(cptr->Q[ivel+dim]) + ss)/gptr->dx[dim];
		if(inv_t == 0)
			inv_t = 1.0/DTMAX;
		tmin = std::min(tmin, 0.5/inv_t);
	}
	return std::min(tmin, DTMAX);
}

void HydroDynamics::calcFluxes() const {
	if(!SWEEPX){
		for(int dim = 0; dim < gptr->ND; ++dim){
			for(GridJoin* jptr = gptr->fjoin[dim]; jptr != NULL; jptr = jptr->next){
				if(gptr->ORDER_S > 0)
					HLLC(jptr->lcell->QR[dim], jptr->F, jptr->rcell->QL[dim], dim);
				else
					HLLC(jptr->lcell->Q, jptr->F, jptr->rcell->Q, dim);
			}
		}
		SWEEPX = !SWEEPX;
	}
	else{
		for(int dim = gptr->ND-1; dim >= 0; --dim){
			for(GridJoin* jptr = gptr->fjoin[dim]; jptr != NULL; jptr = jptr->next){
				if(gptr->ORDER_S > 0)
					HLLC(jptr->lcell->QR[dim], jptr->F, jptr->rcell->QL[dim], dim);
				else
					HLLC(jptr->lcell->Q, jptr->F, jptr->rcell->Q, dim);
			}
		}
		SWEEPX = !SWEEPX;
	}
}
void HydroDynamics::advSolution(const double& dt) const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
		//cptr->U[ihii] = cptr->Q[ihii]*cptr->Q[iden];
		for (int dim = 0; dim < gptr->ND; ++dim) {
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
		//QfromU(cptr->Q, cptr->U);
	}
}

void HydroDynamics::updateBoundaries() const {
	for(int i = 0; i < (int)gptr->boundaries.size(); ++i)
		gptr->boundaries[i]->applyBC();
}

void HydroDynamics::HLLC(double Q_l[], double F[], double Q_r[], const int& dim) const {
	double U_cl[NU], U_cr[NU], U_l[NU], U_r[NU], F_l[NU], F_r[NU], F_cl[NU], F_cr[NU];
	for(int i = 0; i < NU; ++i)
		U_cl[i] = U_cr[i] = U_l[i] = U_r[i] = F_l[i] = F_r[i] = F_cl[i] = F_cr[i] = 0;
	UfromQ(U_l, Q_l);
	UfromQ(U_r, Q_r);
	FfromU(F_l, U_l, dim);
	FfromU(F_r, U_r, dim);

	double S_l = Q_l[ivel+dim] - soundSpeed(Q_l[ipre], Q_l[iden]);
	double S_r = Q_r[ivel+dim] + soundSpeed(Q_r[ipre], Q_r[iden]);

	// S_c = (P_r-P_l+D_l*u_l*(S_l-u_l) - D_r*u_r*(S_r-u_r))/(D_l*(S_l-u_l) - D_r*(S_r-u_r));
	double S_c = (Q_r[ipre]-Q_l[ipre]+Q_l[iden]*Q_l[ivel+dim]*(S_l-Q_l[ivel+dim])-Q_r[iden]*Q_r[ivel+dim]*(S_r-Q_r[ivel+dim]))
			/ (Q_l[iden]*(S_l-Q_l[ivel+dim])-Q_r[iden]*(S_r-Q_r[ivel+dim]));
	double A_r = Q_r[iden]*(S_r-Q_r[ivel+dim])/(S_r-S_c);
	U_cr[iden] = A_r;
	for(int id = 0; id < gptr->ND; ++id)
		U_cr[ivel+id] = A_r*Q_r[ivel+id];
	U_cr[ivel+dim] = A_r*S_c;
	U_cr[ipre] = A_r*((U_r[ipre]/Q_r[iden]) + (S_c-Q_r[ivel+dim])*(S_c + Q_r[ipre]/(Q_r[iden]*(S_r-Q_r[ivel+dim]))));
	U_cr[ihii] = A_r*Q_r[ihii];

	double A_l = Q_l[iden]*(S_l-Q_l[ivel+dim])/(S_l-S_c);
	U_cl[iden] = A_l;
	for(int id = 0; id < gptr->ND; ++id)
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

void HydroDynamics::applySrcTerms(const double& dt, const Radiation* rad) const {
	/* Geometry */
	bool rtp = (gptr->GEOMETRY == SPHERICAL);
	bool rzp = (gptr->GEOMETRY== CYLINDRICAL);
	//bool xyz = (gptr->GEOMETRY == CARTESIAN);
	for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next) {
		/* Geometric source term (if any) */
		if (rzp) {
			double r = gptr->dx[0]*(cptr->xc[0] + 0.5);
			cptr->U[ivel+0] += dt*cptr->Q[ipre]/r;
		}
		else if (rtp){
			double r, area1, area2;
			area1 = cptr->rjoin[0]->area;
			area2 = cptr->ljoin[0]->area;
			r = cptr->vol/(area1 - area2);
			cptr->U[ivel+0] += dt*cptr->Q[ipre]/r;
		}
		//====================//
		/* Radiative Transfer */
		if(rad->nstars > 0){
			//double HII_old, mu_old, T_old;
			double HII_new, mu_new, E_old, E_new, T_new;
			//HII_old = cptr->U[ihii]/cptr->U[iden];
			HII_new = cptr->Q[ihii];
			//mu_old = 1.0/(HII_old + 1.0);
			mu_new = 1.0/(HII_new + 1.0);
			//T_old = (TMIN + (TMAX-TMIN)*HII_old) / (P_SCALE/RHO_SCALE);
			double ke = 0.0;
			for (int dim = 0; dim < gptr->ND; ++dim)
		 	   ke += 0.5*cptr->U[ivel+dim]*cptr->U[ivel+dim]/cptr->U[iden];
			E_old = cptr->U[ipre];
			T_new = (rad->TMIN + (rad->TMAX-rad->TMIN)*HII_new);
			E_new = GAS_CONST*cptr->U[iden]*T_new/(mu_new*(GAMMA - 1.0)) + ke;
			//if(t_hgs > 5*TREC / T_SCALE && rtcoupling)
			if(RTCOUPLING)
				cptr->U[ipre] += (E_new - E_old);
			cptr->U[ihii] = cptr->Q[ihii]*cptr->Q[iden];
			//so[iqtau] += (cptr->pa[iqtau] - cptr->qa[iqtau])/dt;
			//so[iqdtau] += (cptr->pa[iqdtau] - cptr->qa[iqdtau])/dt;
		}
		//====================//
	}
}
void HydroDynamics::fixSolution() const {
	for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next) {
		cptr->U[iden] = std::max(cptr->U[iden], DFLOOR);
		double ke = 0.0;
		for(int dim = 0; dim < gptr->ND; ++dim)
			ke += 0.5*cptr->U[ivel+dim]*cptr->U[ivel+dim]/cptr->U[iden];
		cptr->U[ipre] = std::max(cptr->U[ipre], PFLOOR/(GAMMA - 1.0) + ke);
		cptr->U[ihii] = std::max(std::min(cptr->U[ihii], cptr->U[iden]), 0.0);
	}
}
void HydroDynamics::Qisnan(const int& id, const int& i, const int& xc, const int& yc, const int& zc) const {
	GridCell* cptr = gptr->locate(xc, yc, zc);
	if(std::isnan(cptr->Q[i]))
		std::cerr << id << '\n';
}
