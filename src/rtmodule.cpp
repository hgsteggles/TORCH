/* rtmodule.c */

#include "rtmodule.hpp"
#include "grid3d.hpp"
#include "gridcell.hpp"
#include "partition.hpp"
#include "mpihandler2.hpp"
#include "parameters.hpp"
#include "star.hpp"

#include <math.h>
#include <iostream>
#include <fstream>
#include <cstring>

Radiation::Radiation(const RadiationParameters& rp, Grid3D* grid) : gptr(grid), nstars(0), snapToVertex(false), snapMod(0) {
	K1 = rp.K1;
	K2 = rp.K2;
	K3 = rp.K3;
	K4 = rp.K4;
	P_I_CROSS_SECTION = rp.P_I_CROSS_SECTION;
	ALPHA_B = rp.ALPHA_B;
	TAU_0 = rp.TAU_0;
	SOURCE_S = rp.SOURCE_S;
	NHI = rp.NHI;
	TMIN = rp.TMIN;
	TMAX = rp.TMAX;
	SCHEME = rp.SCHEME;
	H_MASS = rp.H_MASS;
}

void Radiation::addStar(Star src, MPIHandler& mpih, bool snapToVertex) {
	this->snapToVertex = snapToVertex;
	this->snapMod = snapToVertex ? 0.5 : 0.0;
	int dec[3] = {0, 0, 0};
	for (int i = 0; i < gptr->ND && snapToVertex; ++i) {
		if (src.x[i] == gptr->lcell->xc[i]+1)
			dec[i] = -1;
	}
	if (src.x[0] >= gptr->fcell->xc[0] && src.x[0] <= gptr->lcell->xc[0]+1) {
		src.fcausal = gptr->locate((int)(src.x[0])+dec[0], (int)(src.x[1])+dec[1], (int)(src.x[2])+dec[2]);
		src.setCore(mpih.getRank());
	}
	else if (src.x[0] < gptr->fcell->xc[0]) {
		src.fcausal = gptr->locate(gptr->fcell->xc[0], (int)(src.x[1])+dec[1], (int)(src.x[2])+dec[2]);
		src.setCore(mpih.getRank()-1);
	}
	else if (src.x[0] > gptr->lcell->xc[0]+1) {
		src.fcausal = gptr->locate(gptr->lcell->xc[0], (int)(src.x[1])+dec[1], (int)(src.x[2])+dec[2]);
		src.setCore(mpih.getRank()+1);
	}
	stars.push_back(src);
	nstars++;
	gptr->buildCausal(src.fcausal);
}

int Radiation::getRayPlane(GridCell* cptr, const int& starID) const {
	int plane = 0;
	double dxcs_max = 0;
	for(int i = 0; i < gptr->ND; ++i) {
		double dxcs = fabs(cptr->xc[i]+snapMod-stars[starID].x[i]);
		if (dxcs > dxcs_max) {
			dxcs_max = dxcs;
			plane = i;
		}
	}
	return plane;
}

void Radiation::updateTauSC(bool average, GridCell* cptr, const int& starID) const {
	double dist2 = 0;
	for (int i = 0; i < gptr->ND; ++i)
		dist2 += (cptr->xc[i] + snapMod - stars[starID].x[i])*(cptr->xc[i] + snapMod - stars[starID].x[i]);
	if(dist2 > 0.95){
		bool test = false;
		int plane = getRayPlane(cptr, starID);
		int irot[3] = {(plane+1)%3, (plane+2)%3, (plane%3)};
		double d[3] = {0.0, 0.0, 0.0};
		for(int i = 0; i < 3; i++) {
			d[i] = cptr->xc[irot[i]]-stars[starID].x[irot[i]];
			if (irot[i] > gptr->ND - 1) d[i] += snapMod;
		}
		int s[3] = {d[0] < -gptr->dx[0]/10.0 ? -1 : 1, d[1] < -gptr->dx[1]/10.0 ? -1 : 1, d[2] < -gptr->dx[2]/10.0 ? -1 : 1};
		int LR[3] = {abs(d[0]) < gptr->dx[0]/10.0 ? 0 : s[0], abs(d[1]) < gptr->dx[1]/10.0 ? 0 : s[1], abs(d[2]) < gptr->dx[2]/10.0 ? 0 : s[2]};

		double ic[3] = {cptr->xc[irot[0]]-0.5*(s[2]*d[0]/d[2]),
				cptr->xc[irot[1]]-0.5*(s[2]*d[1]/d[2]),
				cptr->xc[irot[2]]-0.5*(s[2])};
		GridCell* e[4] = {gptr->traverse3D(irot[0], irot[1], irot[2], 0, 0, -LR[2], cptr),
				gptr->traverse3D(irot[0], irot[1], irot[2], 0, -LR[1], -LR[2], cptr),
				gptr->traverse3D(irot[0], irot[1], irot[2], -LR[0], 0, -LR[2], cptr),
				gptr->traverse3D(irot[0], irot[1], irot[2], -LR[0], -LR[1], -LR[2], cptr)};
		if (e[0] == NULL && e[1] == NULL && e[2] == NULL && e[3] == NULL) {
			e[0] = gptr->traverseOverJoins3D(irot[0], irot[1], irot[2], 0, 0, -LR[2], cptr);
			e[1] = gptr->traverseOverJoins3D(irot[0], irot[1], irot[2], 0, -LR[1], -LR[2], cptr);
			e[2] = gptr->traverseOverJoins3D(irot[0], irot[1], irot[2], -LR[0], 0, -LR[2], cptr);
			e[3] = gptr->traverseOverJoins3D(irot[0], irot[1], irot[2], -LR[0], -LR[1], -LR[2], cptr);
		}
		double delta[2] = {fabs(2.0*ic[0]-2.0*cptr->xc[irot[0]]+s[0]), fabs(2.0*ic[1]-2.0*cptr->xc[irot[1]]+s[1])};
		double w[4];
		w[0] = (d[0] != 0) || (d[1] != 0) || (d[2] != 0) ? delta[0]*delta[1] : 0;
		w[1] = (d[1] != 0) ? delta[0]*(1.0-delta[1]) : 0;
		w[2] = (d[0] != 0) ? (1.0-delta[0])*delta[1] : 0;
		w[3] = ((d[0] != 0) && (d[1] != 0)) ? (1.0-delta[0])*(1.0-delta[1]) : 0;
		double tau[4] = {0.0, 0.0, 0.0, 0.0};
		double w_raga[4];
		for(int i = 0; i < 4; ++i){
			if (e[i] != NULL)
				tau[i] = e[i]->R[average ? itauta : itau]+e[i]->R[average ? idtauta : idtau];
			w_raga[i] = w[i]/std::max(TAU_0, tau[i]);
		}
		double sum_w = w_raga[0]+w_raga[1]+w_raga[2]+w_raga[3];
		double newtau = 0.0;
		for(int i = 0; i < 4; ++i){
			w_raga[i] = w_raga[i]/sum_w;
			newtau += w_raga[i]*tau[i];
		}
		cptr->R[average ? itauta : itau] = newtau;
		if(test) {
			std::cerr << "d[0] = " << d[0] << std::endl;
			std::cerr << "d[1] = " << d[1] << std::endl;
			std::cerr << "d[2] = " << d[2] << std::endl;
			std::cerr << "ic[0] = " << ic[0] << std::endl;
			std::cerr << "ic[1] = " << ic[1] << std::endl;
			std::cerr << "ic[2] = " << ic[2] << std::endl;
			std::cerr << "delta[0] = " << delta[0] << std::endl;
			std::cerr << "delta[1] = " << delta[1] << std::endl;
		}
	}
	else
		cptr->R[average ? itauta : itau] = 0;
}


void Radiation::calculateNearestNeighbours(const int& starID) const {
	for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->nextcausal) {
		int plane = getRayPlane(cptr, starID);
		int irot[3] = {(plane+1)%3, (plane+2)%3, (plane%3)};
		double d[3] = {0.0, 0.0, 0.0};
		for(int i = 0; i < 3; i++) {
			d[i] = cptr->xc[irot[i]]-stars[starID].x[irot[i]];
			if (irot[i] > gptr->ND - 1) d[i] += snapMod;
		}
		int s[3] = {d[0] < -gptr->dx[0]/10.0 ? -1 : 1, d[1] < -gptr->dx[1]/10.0 ? -1 : 1, d[2] < -gptr->dx[2]/10.0 ? -1 : 1};
		int LR[3] = {std::abs(d[0]) < gptr->dx[0]/10.0 ? 0 : s[0], std::abs(d[1]) < gptr->dx[1]/10.0 ? 0 : s[1], std::abs(d[2]) < gptr->dx[2]/10.0 ? 0 : s[2]};
		cptr->NN[0] = gptr->traverse3D(irot[0], irot[1], irot[2], 0, 0, -LR[2], cptr);
		cptr->NN[1] = gptr->traverse3D(irot[0], irot[1], irot[2], 0, -LR[1], -LR[2], cptr);
		cptr->NN[2] = gptr->traverse3D(irot[0], irot[1], irot[2], -LR[0], 0, -LR[2], cptr);
		cptr->NN[3] = gptr->traverse3D(irot[0], irot[1], irot[2], -LR[0], -LR[1], -LR[2], cptr);
		if (cptr->NN[0] == NULL)
			cptr->NN[0] = gptr->traverseOverJoins3D(irot[0], irot[1], irot[2], 0, 0, -LR[2], cptr);
		if (cptr->NN[1] == NULL)
			cptr->NN[1] = gptr->traverseOverJoins3D(irot[0], irot[1], irot[2], 0, -LR[1], -LR[2], cptr);
		if (cptr->NN[2] == NULL)
			cptr->NN[2] = gptr->traverseOverJoins3D(irot[0], irot[1], irot[2], -LR[0], 0, -LR[2], cptr);
		if (cptr->NN[3] == NULL)
			cptr->NN[3] = gptr->traverseOverJoins3D(irot[0], irot[1], irot[2], -LR[0], -LR[1], -LR[2], cptr);
		double ic[3] = {cptr->xc[irot[0]]-0.5*(s[2]*d[0]/d[2]),	cptr->xc[irot[1]]-0.5*(s[2]*d[1]/d[2]),	cptr->xc[irot[2]]-0.5*(s[2])};
		double delta[2] = {std::abs(2.0*ic[0]-2.0*cptr->xc[irot[0]]+s[0]), std::abs(2.0*ic[1]-2.0*cptr->xc[irot[1]]+s[1])};
		cptr->NN_weights[0] = (d[0] != 0) || (d[1] != 0) || (d[2] != 0) ? delta[0]*delta[1] : 0;
		cptr->NN_weights[1] = (d[1] != 0) ? delta[0]*(1.0-delta[1]) : 0;
		cptr->NN_weights[2] = (d[0] != 0) ? (1.0-delta[0])*delta[1] : 0;
		cptr->NN_weights[3] = ((d[0] != 0) && (d[1] != 0)) ? (1.0-delta[0])*(1.0-delta[1]) : 0;
	}
}

void Radiation::updateTauSC2(bool average, GridCell* cptr, const int& starID) const {
	double dist2 = 0;
	for (int i = 0; i < gptr->ND; ++i)
		dist2 += (cptr->xc[i] + snapMod - stars[starID].x[i])*(cptr->xc[i] + snapMod - stars[starID].x[i]);
	if(dist2 > 0.95){
		double tau[4] = {0.0, 0.0, 0.0, 0.0};
		double w_raga[4];
		for(int i = 0; i < 4; ++i){
			if (cptr->NN[i] != NULL)
				tau[i] = cptr->NN[i]->R[average ? itauta : itau]+cptr->NN[i]->R[average ? idtauta : idtau];
			w_raga[i] = cptr->NN_weights[i]/std::max(TAU_0, tau[i]);
		}
		double sum_w = w_raga[0]+w_raga[1]+w_raga[2]+w_raga[3];
		double newtau = 0.0;
		for(int i = 0; i < 4; ++i){
			w_raga[i] = w_raga[i]/sum_w;
			newtau += w_raga[i]*tau[i];
		}
		cptr->R[average ? itauta : itau] = newtau;
	}
	else
		cptr->R[average ? itauta : itau] = 0;
}

double Radiation::temperature(GridCell* cptr) const {
	return (cptr->Q[ipre]/cptr->Q[iden])*((1.0/(cptr->Q[ihii] + 1.0))/GAS_CONST);
}
double Radiation::alphaB(GridCell* cptr) const {
	//double T = temperature(cptr);
	//double aB = 0.0;
	//double aB = ALPHA_B;//*pow(T/10000, -0.7); // cm^3 s^-1
	return ALPHA_B;
}
double Radiation::cellPathLength(GridCell* cptr, const int& starID) const {
	double denom;
	double d[3] = {0.0, 0.0, 0.0};
	for(int i = 0; i < gptr->ND; ++i)
		d[i] = cptr->xc[i]+snapMod-stars[starID].x[i];
	if(d[2] >= d[1] && d[2] >= d[0])
		denom = d[2];
	else if(d[1] > d[2] && d[1] >= d[0])
		denom = d[1];
	else
		denom = d[0];
	if(denom != 0)
		return sqrt(d[0]*d[0]*gptr->dx[0]*gptr->dx[0]+d[1]*d[1]*gptr->dx[1]*gptr->dx[1]+d[2]*d[2]*gptr->dx[2]*gptr->dx[2])/denom;
	else
		return 0.0;
}
double Radiation::shellVolume(GridCell* cptr, const int& starID) const {
	double r_sqrd = 0;
	for(int i = 0; i < gptr->ND; ++i)
		r_sqrd += (cptr->xc[i]+snapMod-stars[starID].x[i])*(cptr->xc[i]+snapMod-stars[starID].x[i])*gptr->dx[i]*gptr->dx[i];
	//r_min = sqrt(r_sqrd) - 0.5*ds;
	//r_max = sqrt(r_sqrd) + 0.5*ds;
	//if(ND == 1)
	//	return dx*dx*dx;
	//else if(ND == 2)
	//	if(r_sqrd <= 0.1)
	//		return PI*(r_max*r_max - r_min*r_min)*dx;
	//	else
	//		return 2.0*PI*sqrt(r_sqrd)*ds*dx;
	//else if(ND == 3)
	if(r_sqrd <= 0.1) {
		double r_min = sqrt(r_sqrd) - 0.5*cptr->ds;
		double r_max = sqrt(r_sqrd) + 0.5*cptr->ds;
		return 4.0*PI*(r_max*r_max*r_max - r_min*r_min*r_min)/3.0;
	}
	else
		return 4.0*PI*r_sqrd*cptr->ds;
}
double Radiation::PIrate(const double& frac, const double& T, const double& delT, GridCell* cptr, const int& starID) const {
	double n_HI = (1.0-frac)*cptr->Q[iden] / H_MASS;
	if(frac >= 1 || n_HI == 0 || cptr->shellVol == 0)
		return 0.0;
	else {
		return SOURCE_S*exp(-T)*(1.0-exp(-delT))/(n_HI*cptr->shellVol);
	}
}
double Radiation::HIIfracDot(const double& A_pi, const double& frac, GridCell* cptr) const {
	//if(COLLISIONS)
		//A_ci = p.CIlookup[T];
	return (1.0-frac)*A_pi - (frac*frac*alphaB(cptr)*(cptr->Q[iden] / H_MASS)); // + (1-x)*(x*n_H*A_ci);
}
void Radiation::update_dtau(GridCell* cptr, const int& starID) const {
	if(!isStar(cptr)){
		double n_H = cptr->Q[iden] / H_MASS;
		cptr->R[idtau] = (1.0-cptr->Q[ihii])*n_H*P_I_CROSS_SECTION*cptr->ds;
		cptr->R[idtauta] = (1.0-cptr->R[ihiita])*n_H*P_I_CROSS_SECTION*cptr->ds;
	}
}
/*
double Radiation::radHeatCool(double dt, GridCell* cptr) const {
	double T, heat, cool, n_H, n_HII_old, n_HII_new;
	T = temperature(cptr);
	n_HII_old = cptr->U[ihii]/(H_MASS_G/M_SCALE);
	n_HII_new = cptr->Q[ihii]*cptr->Q[iden]/(H_MASS_G/M_SCALE);
	n_H = (cptr->U[iden]/(H_MASS_G/M_SCALE)) * (1.0/(L_SCALE*L_SCALE*L_SCALE)); // cm^-3
	//heat = (HV-H_IE)*(n_HII_new - n_HII_old) / dt;
	heat = 2e-26*n_H / (E_SCALE/(T_SCALE*L_SCALE*L_SCALE*L_SCALE));
	cool = 2e-19*n_H*n_H*(exp((-1.184e5)/(T+1e3)) + 1.4e-9*sqrt(T)*exp(-92/T)); // erg s^-1 cm^-3
	cool /= E_SCALE/(T_SCALE*L_SCALE*L_SCALE*L_SCALE);
	return n_H*heat - cool;
}
*/
void Radiation::doric(const double& dt, double& HII, double& HII_avg, const double& A_pi, GridCell* cptr) const {
	double xeq, inv_ti;
	double aB = alphaB(cptr);
	double epsilon = 1.0e-14;
	//double frac_avg_old;
	double HII_old = HII;
	//frac_avg_old = frac_avg;
	double n_HII_avg = HII_avg*cptr->Q[iden] / H_MASS;
	if(n_HII_avg == 0.0)
		xeq = 1.0;
	else
		xeq = A_pi/(A_pi + n_HII_avg*aB);
	inv_ti = A_pi + n_HII_avg*aB;
	HII = xeq + (HII_old-xeq)*exp(-dt*inv_ti);
	HII = std::max(std::min(HII, 1.0), 0.0);
	if(1.0 - HII < epsilon)
		HII = 1.0 - epsilon;
	if(dt*inv_ti < 1.0e-8)
		HII_avg = HII_old;
	else{
		HII_avg = xeq + (HII_old-xeq)*(1.0-exp(-dt*inv_ti))/(dt*inv_ti);
		HII_avg = std::max(std::min(HII_avg, 1.0), 0.0);
	}
	if(1.0 - HII_avg < epsilon)
		HII_avg = 1.0 - epsilon;
}
void Radiation::update_HIIfrac(const double& dt, GridCell* cptr, const int& starID) const {
	int scheme = 0;
	if(!isStar(cptr)){
		if(scheme == 0){
			double convergence2 = 1.0e-3;
			double convergence_frac = 1.0e-5;
			double HII_avg_old, n_H;
			int niter = 0;
			bool converged = false;
			n_H = cptr->Q[iden] / H_MASS;
			double HII = cptr->Q[ihii];
			//HII_avg = cptr->R[ihiita];
			double HII_avg = HII;
			//tau_avg = cptr->R[itau];
			while(!converged){
				niter++;
				HII_avg_old = HII_avg;
				HII = cptr->Q[ihii];
				double dtau_avg = (1.0-HII_avg)*n_H*P_I_CROSS_SECTION*cptr->ds;
				double A_pi = PIrate(HII_avg, cptr->R[itauta], dtau_avg, cptr, starID);
				doric(dt, HII, HII_avg, A_pi, cptr);
				if((fabs((HII_avg-HII_avg_old)/HII_avg) < convergence2 || (HII_avg < convergence_frac)))
					converged = true;
				if(niter > 5000 || HII != HII){
					std::cout << "ERROR: implicit method not converging." << '\n';
					cptr->printInfo();
					exit(EXIT_FAILURE);
				}
			}
			cptr->Q[ihii] = HII;
			cptr->R[ihiita] = HII_avg;
		}
		else if(scheme == 1){
			double HII, HII_dummy, A_pi, tau, dtau;
			HII = cptr->Q[ihii];
			tau = cptr->R[itau];
			dtau = cptr->R[idtau];
			A_pi = PIrate(HII, tau, dtau, cptr, starID);
			//set_HIIfrac(HII+dt*HIIfracDot(A_pi, HII) );
			HII_dummy = HII;
			doric(dt, HII, HII_dummy, A_pi, cptr);
			cptr->Q[ihii] = HII;
		}
	}
	else{
		cptr->Q[ihii] = 1;
	}
}
double Radiation::getTimeStep(const double& dt_dyn) const {
	static bool once = false;
	double dt = dt_dyn, dt1, dt2, dt3, dt4;
	if (nstars > 0) {
		GridCell* srcptr = stars[0].fcausal;
		for(GridCell* cptr = srcptr; cptr != NULL; cptr = cptr->nextcausal){
			dt1 = dt2 = dt3 = dt4 = dt_dyn;
			if(K1 != 0.0) {
				if (!once) {
					dt1 = K1*1.0/(cptr->Q[iden]*alphaB(cptr)/H_MASS);
					once = true;
				}
				else
					dt1 = K1*1.0/(cptr->Q[iden]*cptr->Q[ihii]*alphaB(cptr)/H_MASS);
			}
			if(K2 != 0.0)
				dt2 = dt_dyn;
			if(K3 != 0.0) {
				double A_pi = PIrate(cptr->Q[ihii], cptr->R[itau], cptr->R[idtau], cptr, 0);
				double fracRate = HIIfracDot(A_pi, cptr->Q[ihii], cptr);
				if (fracRate != 0.0)
					dt3 = K3*std::max(0.05, 1.0 - cptr->Q[ihii])/std::abs(fracRate);
			}
			if(K4 != 0.0) {
				double A_pi = PIrate(cptr->Q[ihii], cptr->R[itau], cptr->R[idtau], cptr, 0);
				double fracRate = HIIfracDot(A_pi, cptr->Q[ihii], cptr);
				if (fracRate != 0.0)
					dt4 = K4*1.0/fabs(fracRate); // timestep criterion [Mackey 2012].
			}
			dt = std::min(dt, std::min(dt_dyn, std::min(dt1, std::min(dt2, std::min(dt3, dt4)))));
			if(dt == 0.0){
				std::cerr << dt << '\t' << dt1 << '\t' << dt_dyn << '\t' << cptr->xc[1] << '\n';
				break;
			}
		}
	}
	return dt;
}
void Radiation::transferRadiation(const double& dt, double& IF) const {
	if (nstars > 0) {
		GridCell* srcptr = stars[0].fcausal;
		bool average = true;
		/** LOOPS OVER CELLS IN GRID CAUSALLY.
		 *
		 */
		for(GridCell* cptr = srcptr; cptr != NULL; cptr = cptr->nextcausal){
			updateTauSC2(average==false, cptr, 0);
			updateTauSC2(average==true, cptr, 0);
			update_HIIfrac(dt, cptr, 0);
			update_dtau(cptr, 0);
		}
	}
}

void Radiation::transferRadiation(const double& dt, double& IF, MPIHandler& mpih) const {
	if (nstars > 0) {
		const unsigned int NELEMENTS = gptr->boundaries[0]->ghostcells.size()*gptr->boundaries[0]->ghostcells[0].size()*NR;
		double* msgArray = new double[NELEMENTS];
		if(stars[0].core != mpih.getRank()) {
			Boundary* part = NULL;
			if(stars[0].core < mpih.getRank()) {
				mpih.receive(msgArray, NELEMENTS, mpih.getRank()-1, RADIATION_MSG);
				part = gptr->boundaries[0];
			}
			else if (stars[0].core > mpih.getRank()) {
				mpih.receive(msgArray, NELEMENTS, mpih.getRank()+1, RADIATION_MSG);
				part = gptr->boundaries[3];
			}
			int id = 0;
			for (int i = 0; i < (int)part->ghostcells.size(); ++i) {
				for (int j = 0; j < (int)part->ghostcells[i].size(); ++j) {
					for(int ir = 0; ir < NR; ++ir)
						part->ghostcells[i][j]->R[ir] = msgArray[id++];
				}
			}
		}
		transferRadiation(dt, IF);
		int face1 = 0, face2 = 3;
		if(mpih.getRank() == 0 || stars[0].core < mpih.getRank())
			face1 = -1;
		if(mpih.getRank() == mpih.nProcessors()-1 || stars[0].core > mpih.getRank())
			face2 = -1;
		for(int i = 0; i < (int)gptr->boundaries.size(); ++i) {
			if(gptr->boundaries[i]->face == face1 || gptr->boundaries[i]->face == face2) {
				Partition* part = (Partition*)gptr->boundaries[i];
				int dim = part->face%3;
				int id = 0;
				for (int i = 0; i < (int)part->ghostcells.size(); ++i) {
					for (int j = 0; j < (int)part->ghostcells[i].size(); j++) {
						GridCell* cptr = NULL;
						if (part->face < 3)
							cptr = part->ghostcells[i][j]->rjoin[dim]->rcell;
						else
							cptr = part->ghostcells[i][j]->ljoin[dim]->lcell;
						for(int ir = 0; ir < NR; ++ir)
							msgArray[id++] = cptr->R[ir];
					}
				}
				mpih.send(msgArray, NELEMENTS, part->destination, RADIATION_MSG);
			}
		}
		delete[] msgArray;
	}
}

void Radiation::rayTrace() const {
	static bool time = 0;
	calculateNearestNeighbours(0);
	if(time == 0 && nstars > 0){
		for(GridCell* cptr = stars[0].fcausal; cptr != NULL; cptr = cptr->nextcausal, ++time){
			cptr->ds = cellPathLength(cptr, 0);
			cptr->shellVol = shellVolume(cptr, 0);
			if(!(!snapToVertex && cptr == stars[0].fcausal))
				update_dtau(cptr, 0);
		}
	}
}

bool Radiation::isStar(GridCell* cptr) const {
	for (int i = 0; i < nstars && !snapToVertex; ++i) {
		if (cptr->xc[0] == stars[i].x[0] && cptr->xc[1] == stars[i].x[1] && cptr->xc[2] == stars[i].x[2])
			return true;
	}
	return false;
}
