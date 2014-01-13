/* hydro.cpp */

#include "hydro.hpp"
#include <typeinfo>

HydroDynamics::HydroDynamics(const HydroParameters& hp) {
	GAMMA = hp.GAMMA;
	DFLOOR = hp.DFLOOR;
	PFLOOR = hp.PFLOOR;
	DTMAX = hp.DTMAX;
}

void HydroDynamics::globalWfromU(Grid3D* gptr) const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->nextCell3D(cptr)){
		for(int iu = 0; iu < NU; iu++)
			cptr->W[iu] = cptr->U[iu];
	}
}
void HydroDynamics::globalUfromW(Grid3D* gptr) const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->nextCell3D(cptr)){
		for(int iu = 0; iu < NU; iu++)
			cptr->U[iu] = cptr->W[iu];
	}
}
void HydroDynamics::globalQfromU(Grid3D* gptr) const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->nextCell3D(cptr))
		QfromU(cptr->Q, cptr->U);
}
void HydroDynamics::globalUfromQ(Grid3D* gptr) const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->nextCell3D(cptr))
		UfromQ(cptr->U, cptr->Q);
}
void HydroDynamics::reconstruct(Grid3D* gptr) const {
	double dQdr1[NU], dQdr2[NU], dQdr[NU];
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->nextCell3D(cptr)){
		for(int dim = 0; dim < ND; dim++){
			for(int iq = 0; iq < NU; iq++){
				if(cptr->left[dim] != NULL)
					dQdr1[iq] = (cptr->Q[iq] - cptr->left[dim]->Q[iq])/gptr->dx[dim];
				else if(cptr->left[dim] == NULL){
					dQdr1[iq] = cptr->Q[iq];
					dQdr1[iq] = cptr->bd[dim]->ghosts[0]->Q[iq];
					dQdr1[iq] = gptr->dx[dim];
					dQdr1[iq] = (cptr->Q[iq] - cptr->bd[dim]->ghosts[0]->Q[iq])/gptr->dx[dim];
				}
				if(cptr->right[dim] != NULL)
					dQdr2[iq] = (cptr->right[dim]->Q[iq] - cptr->Q[iq])/gptr->dx[dim];
				else if(cptr->right[dim] == NULL)
					dQdr2[iq] = (cptr->bd[dim]->ghosts[0]->Q[iq] - cptr->Q[iq])/gptr->dx[dim];
				dQdr[iq] = av(dQdr1[iq], dQdr2[iq]);
				if(iq == ihii && isnan(dQdr[iq])) {
					cerr << "dQdr1 = " << dQdr1[iq] << endl;
					cerr << "dQdr2 = " << dQdr2[iq] << endl;
				}
				cptr->QL[dim][iq] = cptr->Q[iq] - 0.5*(gptr->dx[dim])*dQdr[iq];
				cptr->QR[dim][iq] = cptr->Q[iq] + 0.5*(gptr->dx[dim])*dQdr[iq];
			}
			UfromQ(cptr->UL[dim], cptr->QL[dim]);
			UfromQ(cptr->UR[dim], cptr->QR[dim]);
		}
	}
	for(Boundary* bptr = gptr->fboundary; bptr != NULL; bptr = bptr->next){
		int dim = bptr->face%3;
		for(int iq = 0; iq < NU; iq++){
			dQdr1[iq] = (bptr->ghosts[0]->Q[iq] - bptr->ghosts[1]->Q[iq])/gptr->dx[dim];
			dQdr2[iq] = (bptr->gridcell->Q[iq] - bptr->ghosts[0]->Q[iq])/gptr->dx[dim];
			if(bptr->face > 2){
				dQdr1[iq] = -dQdr1[iq];
				dQdr2[iq] = -dQdr2[iq];
			}
			dQdr[iq] = av(dQdr1[iq], dQdr2[iq]);
			bptr->ghosts[0]->QL[dim][iq] = bptr->ghosts[0]->Q[iq] - 0.5*(gptr->dx[dim])*dQdr[iq];
			bptr->ghosts[0]->QR[dim][iq] = bptr->ghosts[0]->Q[iq] + 0.5*(gptr->dx[dim])*dQdr[iq];
		}
		for(int i = 0; i < 3; i++){
			UfromQ(bptr->ghosts[0]->UL[i], bptr->ghosts[0]->QL[i]);
			UfromQ(bptr->ghosts[0]->UR[i], bptr->ghosts[0]->QR[i]);
		}
	}
}
double HydroDynamics::av(double a, double b) const {
	double result = 0;
	double nom = a*a*b + a*b*b;
	double denom = a*a + b*b;
	if(denom > numeric_limits<double>::min())
		result = nom/denom;
	return result;
}

void HydroDynamics::UfromQ(double u[], double q[]) const {
	double ke = 0;
	u[iden] = q[iden] = max(q[iden], DFLOOR);
	for(int dim = 0; dim < 3; dim++){
		u[ivel+dim] = q[ivel+dim]*q[iden];
		ke += 0.5*q[iden]*q[ivel+dim]*q[ivel+dim];
	}
	q[ipre] = max(q[ipre], PFLOOR);
	u[ipre] = q[ipre]/(GAMMA - 1.0) + ke;
	q[ihii] = max(min(q[ihii], 1.0), 0.0);
	u[ihii] = q[ihii]*q[iden];
}
void HydroDynamics::QfromU(double q[], double u[]) const {
	q[iden] = u[iden] = max(u[iden], DFLOOR);
	double ke = 0;
	for(int dim = 0; dim < 3; dim++){
		q[ivel+dim] = u[ivel+dim]/u[iden];
		ke += 0.5*u[ivel+dim]*u[ivel+dim]/u[iden];
	}
	u[ipre] = max(u[ipre], PFLOOR/(GAMMA - 1.0) + ke);
	q[ipre] = (u[ipre] - ke)*(GAMMA - 1.0);
	u[ihii] = max(min(u[ihii], u[iden]), 0.0);
	q[ihii] = u[ihii]/u[iden];
}
void HydroDynamics::FfromU(double f[], double u[], int dim) const {
	double ke = 0, pressure;
	for(int id = 0; id < 3; id++)
		ke += 0.5*u[ivel+id]*u[ivel+id]/u[iden];
	pressure = (u[ipre] - ke)*(GAMMA - 1.0);
	f[iden] = u[ivel+dim];
	f[ivel+dim] = u[ivel+dim]*u[ivel+dim]/u[iden] + pressure;
	f[ivel+(dim+1)%3] = u[ivel+(dim+1)%3]*u[ivel+dim]/u[iden];
	f[ivel+(dim+2)%3] = u[ivel+(dim+2)%3]*u[ivel+dim]/u[iden];
	f[ipre] = u[ivel+dim]*(u[ipre] + pressure)/u[iden];
	f[ihii] = u[ivel+dim]*u[ihii]/u[iden];
}
double HydroDynamics::soundSpeed(double pre, double den) const {
	return sqrt(GAMMA*pre/den);
}
double HydroDynamics::CFL(Grid3D* gptr) const {
	double tmin = DTMAX;
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->nextCell3D(cptr)){
		double inv_t = 0;
		double ss = soundSpeed(cptr->Q[iden], cptr->Q[ipre]);
		for(int dim = 0; dim < 3; dim++)
			inv_t += (fabs(cptr->Q[ivel+dim]) + ss)/gptr->dx[dim];
		if(inv_t == 0)
			inv_t = 1.0/DTMAX;
		tmin = min(tmin, 0.5/inv_t);
	}
	return min(tmin, DTMAX);
}

void HydroDynamics::calcFluxes(Grid3D* gptr) const {
	if(!SWEEPX){
		for(int dim = 0; dim < ND; dim++){
			for(GridJoin* jptr = gptr->fjoin[dim]; jptr != NULL; jptr = jptr->next){
				if(ORDER > 0)
					HLLC(jptr->lcell->UR[dim], jptr->F, jptr->rcell->UL[dim], dim);
				else
					HLLC(jptr->lcell->U, jptr->F, jptr->rcell->U, dim);
			}
		}
		SWEEPX = !SWEEPX;
	}
	else{
		for(int dim = ND-1; dim >= 0; dim--){
			for(GridJoin* jptr = gptr->fjoin[dim]; jptr != NULL; jptr = jptr->next){
				if(ORDER > 0)
					HLLC(jptr->lcell->UR[dim], jptr->F, jptr->rcell->UL[dim], dim);
				else
					HLLC(jptr->lcell->U, jptr->F, jptr->rcell->U, dim);
			}
		}
		SWEEPX = !SWEEPX;
	}
}
void HydroDynamics::calcBoundaryFluxes(Grid3D* gptr) const {
	for(Boundary* bptr = gptr->fboundary; bptr != NULL; bptr = bptr->next){
		int dim = bptr->face%3;
		if(bptr->face > 2){
			if(ORDER > 0)
				HLLC(bptr->gridcell->UR[dim], bptr->F, bptr->ghosts[0]->UL[dim], dim);
			else
				HLLC(bptr->gridcell->U, bptr->F, bptr->ghosts[0]->U, dim);
		}
		else{
			if(ORDER > 0)
				HLLC(bptr->ghosts[0]->UR[dim], bptr->F, bptr->gridcell->UL[dim], dim);
			else
				HLLC(bptr->ghosts[0]->U, bptr->F, bptr->gridcell->U, dim);
		}
	}
}
void HydroDynamics::advSolution(double dt, Grid3D* gptr) const {
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->nextCell3D(cptr)){
		//cptr->U[ihii] = cptr->Q[ihii]*cptr->Q[iden];
		for (int dim = 0; dim < ND; dim++) {
			for (int i = 0; i < NU; i++) {
				if(cptr->ljoin[dim] != NULL && cptr->rjoin[dim] != NULL) {
					if(isnan(cptr->ljoin[dim]->area))
						cerr << "here1" << endl;
					if(isnan(cptr->rjoin[dim]->area))
						cerr << "here2" << endl;
					if(isnan(cptr->vol))
						cerr << "here3" << endl;
					if(isnan(cptr->ljoin[dim]->F[i]))
						cerr << "here4" << endl;
					if(isnan(cptr->rjoin[dim]->F[i]))
						cerr << "here5" << endl;
					cptr->U[i] += (dt*cptr->ljoin[dim]->area/cptr->vol)*cptr->ljoin[dim]->F[i];
					cptr->U[i] -= (dt*cptr->rjoin[dim]->area/cptr->vol)*cptr->rjoin[dim]->F[i];
				}
				else if(cptr->ljoin[dim] == NULL && cptr->rjoin[dim] != NULL) {
					cptr->U[i] += (dt*cptr->bd[dim]->area/cptr->vol)*cptr->bd[dim]->F[i];
					cptr->U[i] -= (dt*cptr->rjoin[dim]->area/cptr->vol)*cptr->rjoin[dim]->F[i];
				}
				else if(cptr->ljoin[dim] != NULL && cptr->rjoin[dim] == NULL) {
					cptr->U[i] += (dt*cptr->ljoin[dim]->area/cptr->vol)*cptr->ljoin[dim]->F[i];
					cptr->U[i] -= (dt*cptr->bd[dim]->area/cptr->vol)*cptr->bd[dim]->F[i];
				}
			}
		}
		QfromU(cptr->Q, cptr->U);
	}
}

void HydroDynamics::updateBoundaries(Grid3D* gptr) const {
	for(Boundary* bptr = gptr->fboundary; bptr != NULL; bptr = bptr->next)
		bptr->applyBC();
}

void HydroDynamics::HLLC(double U_l[], double F[], double U_r[], int dim) const {
	double a_l = 0, a_r = 0, S_l = 0, S_c = 0, S_r = 0, A_l = 0, A_r = 0;
	double U_cl[NU], U_cr[NU], Q_l[NU], Q_r[NU];
	double F_l[NU], F_r[NU], F_cl[NU], F_cr[NU];
	for(int i = 0; i < NU; i++)
		U_cl[i] = U_cr[i] = Q_l[i] = Q_r[i] = F_l[i] = F_r[i] = F_cl[i] = F_cr[i] = 0;
	QfromU(Q_l, U_l);
	QfromU(Q_r, U_r);
	FfromU(F_l, U_l, dim);
	FfromU(F_r, U_r, dim);

	a_l = soundSpeed(Q_l[ipre], Q_l[iden]);
	a_r = soundSpeed(Q_r[ipre], Q_r[iden]);
	S_l = Q_l[ivel+dim] - a_l;
	S_r = Q_r[ivel+dim] + a_r;

	double nom1 = Q_r[ipre]-Q_l[ipre]+Q_l[iden]*Q_l[ivel+dim]*(S_l-Q_l[ivel+dim]);
	double nom2 = -Q_r[iden]*Q_r[ivel+dim]*(S_r-Q_r[ivel+dim]);
	double denom = Q_l[iden]*(S_l-Q_l[ivel+dim])-Q_r[iden]*(S_r-Q_r[ivel+dim]);
	S_c = (nom1+nom2)/denom;
	A_r = Q_r[iden]*(S_r-Q_r[ivel+dim])/(S_r-S_c);
	U_cr[iden] = A_r;
	for(int id = 0; id < ND; id++)
		U_cr[ivel+id] = A_r*Q_r[ivel+id];
	U_cr[ivel+dim] = A_r*S_c;
	U_cr[ipre] = A_r*((U_r[ipre]/Q_r[iden]) + (S_c-Q_r[ivel+dim])*(S_c + Q_r[ipre]/(Q_r[iden]*(S_r-Q_r[ivel+dim]))));
	U_cr[ihii] = A_r*Q_r[ihii];

	A_l = Q_l[iden]*(S_l-Q_l[ivel+dim])/(S_l-S_c);
	U_cl[iden] = A_l;
	for(int id = 0; id < ND; id++)
		U_cl[ivel+id] = A_l*Q_l[ivel+id];
	U_cl[ivel+dim] = A_l*S_c;
	U_cl[ipre] = A_l*((U_l[ipre]/Q_l[iden]) + (S_c-Q_l[ivel+dim])*(S_c+Q_l[ipre]/(Q_l[iden]*(S_l-Q_l[ivel+dim]))));
	U_cl[ihii] = A_l*Q_l[ihii];
	for(int i = 0; i < NU; i++){
		F_cl[i] = F_l[i]+S_l*(U_cl[i]-U_l[i]);
		F_cr[i] = F_r[i]+S_r*(U_cr[i]-U_r[i]);
	}
	if(S_l >= 0){
		for(int i = 0; i < NU; i++)
			F[i] = F_l[i];
	}
	else if(S_l <= 0 && S_c >= 0){
		for(int i = 0; i < NU; i++)
			F[i] = F_cl[i];
	}
	else if (S_c <= 0 && S_r >= 0){
		for(int i = 0; i < NU; i++)
			F[i] = F_cr[i];
	}
	else if(S_r <= 0){
		for(int i = 0; i < NU; i++)
			F[i] = F_r[i];
	}
	if(isnan(F[ihii])) {
		for(int i = 0; i < NU; i++){
			cerr << "Q_l[" << i << "] = " << Q_l[i] << endl;
			cerr << "Q_r[" << i << "] = " << Q_r[i] << endl;
		}
		for(int i = 0; i < NU; i++){
			cerr << "U_l[" << i << "] = " << U_l[i] << endl;
			cerr << "U_r[" << i << "] = " << U_r[i] << endl;
		}
		for(int i = 0; i < NU; i++){
			cerr << "F_l[" << i << "] = " << F_l[i] << endl;
			cerr << "F_r[" << i << "] = " << F_r[i] << endl;
		}
		for(int i = 0; i < NU; i++){
			cerr << "U_cl[" << i << "] = " << U_cl[i] << endl;
			cerr << "U_cr[" << i << "] = " << U_cr[i] << endl;
		}
		cerr << "a_l = " << a_l << endl;
		cerr << "a_r = " << a_r << endl;
		cerr << "S_l = " << S_l << endl;
		cerr << "S_r = " << S_r << endl;
		cerr << "S_c = " << S_c << endl;
		for(int iu = 0; iu < NU; iu++)
			cerr << "F[" << iu << "] = " << F[iu] << endl;
		exit(EXIT_FAILURE);
	}
}

void HydroDynamics::applySrcTerms(double dt, Grid3D* gptr, const Radiation& rad) const {
	/* Geometry */
	bool rtp = (gptr->GEOMETRY == SPHERICAL);
	bool rzp = (gptr->GEOMETRY== CYLINDRICAL);
	//bool xyz = (gptr->GEOMETRY == CARTESIAN);
	for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->nextCell3D(cptr)) {
		/* Geometric source term (if any) */
		if (rzp) {
			double r = gptr->dx[0]*(cptr->xc[0] - 0.5);
			cptr->U[ivel+0] += dt*cptr->Q[ipre]/r;
		}
		else if (rtp){
			double r, area1, area2;
			if(cptr->rjoin[0] == NULL)
				area1 = cptr->bd[0]->area;
			else
				area1 = cptr->rjoin[0]->area;
			if(cptr->ljoin[0] == NULL)
				area2 = cptr->bd[0]->area;
			else
				area2 = cptr->ljoin[0]->area;
			r = cptr->vol/(area1 - area2);
			cptr->U[ivel+0] += dt*cptr->Q[ipre]/r;
		}
		//====================//
		/* Radiative Transfer */
		if(gptr->fsrc != NULL){
			//double HII_old, mu_old, T_old;
			double HII_new, mu_new, E_old, E_new, T_new;
			//HII_old = cptr->U[ihii]/cptr->U[iden];
			HII_new = cptr->Q[ihii];
			//mu_old = 1.0/(HII_old + 1.0);
			mu_new = 1.0/(HII_new + 1.0);
			//T_old = (TMIN + (TMAX-TMIN)*HII_old) / (P_SCALE/RHO_SCALE);
			double ke = 0.0;
			for (int dim = 0; dim < ND; dim++)
		 	   ke += 0.5*cptr->U[ivel+dim]*cptr->U[ivel+dim]/cptr->U[iden];
			E_old = cptr->U[ipre];
			T_new = (rad.TMIN + (rad.TMAX-rad.TMIN)*HII_new);
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
void HydroDynamics::fixSolution(Grid3D* gptr) const {
	for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->nextCell3D(cptr)) {
		cptr->U[iden] = max(cptr->U[iden], DFLOOR);
		double ke = 0.0;
		for(int dim = 0; dim < 3; dim++)
			ke += 0.5*cptr->U[ivel+dim]*cptr->U[ivel+dim]/cptr->U[iden];
		cptr->U[ipre] = max(cptr->U[ipre], PFLOOR/(GAMMA - 1.0) + ke);
		cptr->U[ihii] = max(min(cptr->U[ihii], cptr->U[iden]), 0.0);
	}
}
void HydroDynamics::Qisnan(int id, int i, int xc, int yc, int zc, Grid3D* gptr) const {
	GridCell* cptr = gptr->locate(xc, yc, zc);
	if(isnan(cptr->Q[i]))
		cerr << id << endl;
}
