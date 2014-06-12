/**
 * @file radiation.cpp
 */

#include "radiation.hpp"
#include "grid3d.hpp"
#include "gridcell.hpp"
#include "partition.hpp"
#include "mpihandler.hpp"
#include "parameters.hpp"
#include "star.hpp"
#include "thermodynamics.hpp"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>

Radiation::Radiation(const RadiationParameters& rp, const Scalings& sc, Grid3D* grid) : gptr(grid), nstars(0) {
	rparams = new RadiationParameters(rp);
	scale = new Scalings(sc);
}

int Radiation::getRayPlane(GridCell* cptr, const int& starID) const {
	double mod[3] = {stars[starID].mod[0], stars[starID].mod[1], stars[starID].mod[2]};
	int plane = 0;
	double dxcs_max = 0;
	for(int i = 0; i < gptr->gparams->ND; ++i) {
		double dxcs = fabs(cptr->xc[i]+mod[i]-stars[starID].x[i]);
		if (dxcs > dxcs_max) {
			dxcs_max = dxcs;
			plane = i;
		}
	}
	return plane;
}

void Radiation::calculateNearestNeighbours(const int& starID) const {
	double mod[3] = {stars[starID].mod[0], stars[starID].mod[1], stars[starID].mod[2]};
	for (GridCell* cptr = stars[starID].fcausal; cptr != NULL; cptr = cptr->nextcausal) {
		int plane = getRayPlane(cptr, starID);
		int irot[3] = {(plane+1)%3, (plane+2)%3, (plane%3)};
		double d[3] = {0.0, 0.0, 0.0};
		for(int i = 0; i < 3; i++)
			d[i] = cptr->xc[irot[i]] + mod[irot[i]] -stars[starID].x[irot[i]];
		int s[3] = {d[0] < -1.0/10.0 ? -1 : 1, d[1] < -1.0/10.0 ? -1 : 1, d[2] < -1.0/10.0 ? -1 : 1};
		int LR[3] = {std::abs(d[0]) < 1.0/10.0 ? 0 : s[0], std::abs(d[1]) < 1.0/10.0 ? 0 : s[1], std::abs(d[2]) < 1.0/10.0 ? 0 : s[2]};
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
		cptr->NN_weights[0] = (std::abs(d[2]) > 0.9) ? delta[0]*delta[1] : 0;
		cptr->NN_weights[1] = ((std::abs(d[1]) > 0.9) && (std::abs(d[2]) > 0.9)) ? delta[0]*(1.0-delta[1]) : 0;
		cptr->NN_weights[2] = ((std::abs(d[0]) > 0.9) && (std::abs(d[2]) > 0.9)) ? (1.0-delta[0])*delta[1] : 0;
		cptr->NN_weights[3] = ((std::abs(d[0]) > 0.9) && (std::abs(d[1]) > 0.9) && (std::abs(d[2]) > 0.9)) ? (1.0-delta[0])*(1.0-delta[1]) : 0;
	}
}

void Radiation::rayTrace() const {
	static bool time = 0;
	calculateNearestNeighbours(0);
	if(time == 0 && nstars > 0){
		for(GridCell* cptr = stars[0].fcausal; cptr != NULL; cptr = cptr->nextcausal, ++time){
		//for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
			cptr->ds = cellPathLength(cptr, 0);
			cptr->shellVol = shellVolume(cptr, 0);
			update_dtau(cptr, 0);
		}
	}
}

double Radiation::cellPathLength(GridCell* cptr, const int& starID) const {
	double denom;
	double d[3] = {0.0, 0.0, 0.0};
	for(int i = 0; i < gptr->gparams->ND; ++i)
		d[i] = std::abs(cptr->xc[i] + stars[starID].mod[i] - stars[starID].x[i]);
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
	for(int i = 0; i < gptr->gparams->ND; ++i)
		r_sqrd += (cptr->xc[i]+stars[starID].mod[i]-stars[starID].x[i])*(cptr->xc[i]+stars[starID].mod[i]-stars[starID].x[i])*gptr->dx[i]*gptr->dx[i];
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

void Radiation::doric(const double& dt, double& HII, double& HII_avg, const double& A_pi, GridCell* cptr) const {
	double xeq, inv_ti;
	double aB = alphaB(temperature(cptr));
	double epsilon = 1.0e-8;
	//double frac_avg_old;
	double HII_old = HII;
	//frac_avg_old = frac_avg;
	double n_HII_avg = HII_avg*cptr->Q[iden] / rparams->H_MASS;
	if(n_HII_avg == 0.0)
		xeq = 1.0;
	else
		xeq = A_pi/(A_pi + n_HII_avg*aB);
	inv_ti = A_pi + n_HII_avg*aB;
	if (dt*inv_ti == 0) {
	//if(dt*inv_ti < 1.0e-12) {
		/*
		if (inv_ti > 1.0e-7) {
			double HII_exact = xeq + (HII_old-xeq)*exp(-dt*inv_ti);
			double HII_avg_exact = HII_avg = xeq + (HII_old-xeq)*(1.0-exp(-dt*inv_ti))/(dt*inv_ti);
			double HII_error = std::abs(HII_old - HII_exact);
			double HII_avg_error = std::abs(HII_old - HII_avg_exact);
			std::cerr << HII_error << " " << HII_avg_error << '\n';
		}
		*/
		HII = HII_old;
		HII = std::max(std::min(HII, 1.0), 0.0);
		if(1.0 - HII < epsilon)
			HII = 1.0 - epsilon;
		HII_avg = HII;
	}
	else {
		HII = xeq + (HII_old-xeq)*exp(-dt*inv_ti);
		HII = std::max(std::min(HII, 1.0), 0.0);
		if(1.0 - HII < epsilon)
			HII = 1.0 - epsilon;
		HII_avg = xeq + (HII_old-xeq)*(1.0-exp(-dt*inv_ti))/(dt*inv_ti);
		HII_avg = std::max(std::min(HII_avg, 1.0), 0.0);
		if(1.0 - HII_avg < epsilon)
			HII_avg = 1.0 - epsilon;
	}
}

double Radiation::temperature(GridCell* cptr) const {
	return (cptr->Q[ipre]/cptr->Q[iden])*(scale->P/scale->RHO)*((1.0/(cptr->Q[ihii] + 1.0))/GAS_CONST);
}

double Radiation::alphaB(const double& T) const {
	//double T = temperature(cptr);
	//double aB = 0.0;
	//double aB = ALPHA_B;//*pow(T/10000, -0.7); // cm^3 s^-1
	return rparams->ALPHA_B;
}

double Radiation::collIonRate(const double& Temperature) {
	//Fitting formulae from Voronov (1997) ADANDT, 65, 1. (Atomic Data And Nuclear Data Tables)
	double T = Temperature / (scale->P/scale->RHO);
	if (T < 5.0e3) return 0.0;
	double A = 2.91e-8, X = 0.232, K = 0.39;
	int PP = 0;
	double U = 13.6/BOLTZMANN*ERG2EV/T;
	return A*(1.+PP*sqrt(U))*exp(K*log(U) -U)/(X+U); //cm^3/s
}

double Radiation::PIrate(const double& frac, const double& T, const double& delT, GridCell* cptr, const int& starID) const {
	double n_HI = (1.0-frac)*cptr->Q[iden] / rparams->H_MASS;
	if(frac >= 1 || n_HI == 0 || cptr->shellVol == 0)
		return 0.0;
	else {
		return rparams->SOURCE_S*exp(-T)*(1.0-exp(-delT))/(n_HI*cptr->shellVol);
	}
}

double Radiation::HIIfracDot(const double& A_pi, const double& frac, GridCell* cptr) const {
	//if(COLLISIONS)
		//A_ci = p.CIlookup[T];
	return (1.0-frac)*A_pi - (frac*frac*alphaB(temperature(cptr))*(cptr->Q[iden] / rparams->H_MASS)); // + (1-x)*(x*n_H*A_ci);
}

double Radiation::calcTimeStep(const double& dt_dyn) const {
	static bool once = false;
	double dt = dt_dyn, dt1, dt2, dt3, dt4;
	if (nstars > 0) {
		GridCell* srcptr = stars[0].fcausal;
		for(GridCell* cptr = srcptr; cptr != NULL; cptr = cptr->nextcausal){
			dt1 = dt2 = dt3 = dt4 = dt_dyn;
			if(rparams->K1 != 0.0) {
				if (!once) {
					dt1 = rparams->K1*1.0/(cptr->Q[iden]*alphaB(temperature(cptr))/rparams->H_MASS);
					once = true;
				}
				else
					dt1 = rparams->K1*1.0/(cptr->Q[iden]*cptr->Q[ihii]*alphaB(temperature(cptr))/rparams->H_MASS);
				if (dt1 < 1.0e-6) {
					std::cout << "ERROR: rad step too small. Density = " << cptr->Q[iden];
					std::cout << ", Hii = " << cptr->Q[ihii] << '\n';
				}
			}
			if(rparams->K2 != 0.0)
				dt2 = dt_dyn;
			if(rparams->K3 != 0.0) {
				double A_pi = PIrate(cptr->Q[ihii], cptr->R[itau], cptr->R[idtau], cptr, 0);
				double fracRate = HIIfracDot(A_pi, cptr->Q[ihii], cptr);
				if (fracRate != 0.0)
					dt3 = rparams->K3*std::max(0.05, 1.0 - cptr->Q[ihii])/std::fabs(fracRate);
			}
			if(rparams->K4 != 0.0) {
				double A_pi = PIrate(cptr->Q[ihii], cptr->R[itau], cptr->R[idtau], cptr, 0);
				double fracRate = HIIfracDot(A_pi, cptr->Q[ihii], cptr);
				if (fracRate != 0.0)
					dt4 = rparams->K4*1.0/fabs(fracRate); // timestep criterion [Mackey 2012].
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

void Radiation::updateTauSC(bool average, GridCell* cptr, const int& starID) const {
	double mod[3] = {stars[starID].mod[0], stars[starID].mod[1], stars[starID].mod[2]};
	double dist2 = 0;
	for (int i = 0; i < gptr->gparams->ND; ++i)
		dist2 += (cptr->xc[i] + mod[i] - stars[starID].x[i])*(cptr->xc[i] + mod[i] - stars[starID].x[i]);
	if(dist2 > 0.95){
		bool test = false;
		int plane = getRayPlane(cptr, starID);
		int irot[3] = {(plane+1)%3, (plane+2)%3, (plane%3)};
		double d[3] = {0.0, 0.0, 0.0};
		for(int i = 0; i < 3; i++)
			d[i] = cptr->xc[irot[i]] + mod[irot[i]] -stars[starID].x[irot[i]];
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
			w_raga[i] = w[i]/std::max(rparams->TAU_0, tau[i]);
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

void Radiation::updateTauSC2(bool average, GridCell* cptr, const double& dist2) const {
	if(dist2 > 0.95){
		double tau[4] = {0.0, 0.0, 0.0, 0.0};
		double w_raga[4];
		for(int i = 0; i < 4; ++i) {
			if (cptr->NN[i] != NULL)
				tau[i] = cptr->NN[i]->R[average ? itauta : itau]+cptr->NN[i]->R[average ? idtauta : idtau];
			w_raga[i] = cptr->NN_weights[i]/std::max(rparams->TAU_0, tau[i]);
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

void Radiation::updateColDen(GridCell* cptr, const double& dist2) const {
	if(dist2 > 0.95*0.95){
		double colden[4] = {0.0, 0.0, 0.0, 0.0};
		double w_raga[4];
		for(int i = 0; i < 4; ++i) {
			if (cptr->NN[i] != NULL) {
				double dcolden;
				if (dist2 > 1.8*1.8)
					dcolden = (cptr->NN[i]->Q[iden] / rparams->H_MASS)*rparams->P_I_CROSS_SECTION*cptr->NN[i]->ds;
				else if (dist2 > 1.5*1.5) // e.g. cell@{1,1,1}, star@{0,0,0}
					dcolden = (cptr->NN[i]->Q[iden] / rparams->H_MASS)*rparams->P_I_CROSS_SECTION*1.732050808*0.5;
				else if (dist2 > 1.1*1.1) // e.g. cell@{1,1,0}, star@{0,0,0}
					dcolden = (cptr->NN[i]->Q[iden] / rparams->H_MASS)*rparams->P_I_CROSS_SECTION*1.414213562*0.5;
				else // e.g. cell@{1,0,0}, star@{0,0,0}
					dcolden = (cptr->NN[i]->Q[iden] / rparams->H_MASS)*rparams->P_I_CROSS_SECTION*0.5;
				colden[i] = cptr->NN[i]->R[icolden]+dcolden;
			}
			w_raga[i] = cptr->NN_weights[i]/std::max(rparams->TAU_0, colden[i]);
		}
		double sum_w = w_raga[0]+w_raga[1]+w_raga[2]+w_raga[3];
		double newcolden = 0.0;
		for(int i = 0; i < 4; ++i){
			w_raga[i] = w_raga[i]/sum_w;
			newcolden += w_raga[i]*colden[i];
		}
		cptr->R[icolden] = newcolden;
	}
	else
		cptr->R[icolden] = 0;
}

void Radiation::update_dtau(GridCell* cptr, const int& starID) const {
	double n_H = cptr->Q[iden] / rparams->H_MASS;
	cptr->R[idtau] = (1.0-cptr->Q[ihii])*n_H*rparams->P_I_CROSS_SECTION*cptr->ds;
	if (rparams->SCHEME == IMPLICIT)
		cptr->R[idtauta] = (1.0-cptr->R[ihiita])*n_H*rparams->P_I_CROSS_SECTION*cptr->ds;
}

void Radiation::update_HIIfrac(const double& dt, GridCell* cptr, const int& starID) const {
	if(!isStar(cptr)){
		double A_pi = 0;
		double HII = cptr->Q[ihii];
		if(rparams->SCHEME == IMPLICIT || rparams->SCHEME == IMPLICIT2){
			double convergence2 = 1.0e-4;
			double convergence_frac = 1.0e-6;
			double HII_avg_old, n_H;
			int niter = 0;
			bool converged = false;
			n_H = cptr->Q[iden] / rparams->H_MASS;
			//HII_avg = cptr->R[ihiita];
			double HII_avg = HII;
			//tau_avg = cptr->R[itau];
			while(!converged){
				niter++;
				HII_avg_old = HII_avg;
				HII = cptr->Q[ihii];
				double dtau_avg = (1.0-HII_avg)*n_H*rparams->P_I_CROSS_SECTION*cptr->ds;
				A_pi = PIrate(HII_avg, cptr->R[itauta], dtau_avg, cptr, starID);
				doric(dt, HII, HII_avg, A_pi, cptr);
				if((fabs((HII_avg-HII_avg_old)/HII_avg) < convergence2 || (HII_avg < convergence_frac)))
					converged = true;
				if(niter > 5000) {
					std::cout << "HII_avg = " << HII_avg << ", HII_avg_old = " << HII_avg_old << ", HII = " << HII << '\n';
				}
				if(niter > 50000 || HII != HII) {
					std::cout << "ERROR: implicit method not converging." << '\n';
					cptr->printInfo();
					exit(EXIT_FAILURE);
				}
			}
			cptr->R[ihiita] = HII_avg;
		}
		else if(rparams->SCHEME == EXPLICIT){
			double HII_dummy, tau, dtau;
			tau = cptr->R[itau];
			dtau = cptr->R[idtau];
			A_pi = PIrate(HII, tau, dtau, cptr, starID);
			//set_HIIfrac(HII+dt*HIIfracDot(A_pi, HII) );
			HII_dummy = HII;
			doric(dt, HII, HII_dummy, A_pi, cptr);
		}
		//Photoionization heating.
		//double nh = cptr->Q[iden] / rparams->H_MASS;
		//cptr->R[icool] = -(1.0-cptr->Q[ihii])*A_pi*nh*(stars[starID].photEnergy-13.6*EV2ERG/scale->E);
		//Recombination cooling.
		//cptr->R[icool] += calcRecombinationEnergy(temperature(cptr))*cptr->Q[ihii]*cptr->Q[ihii]*nh;
		cptr->R[icool] = 0;
		cptr->Q[ihii] = HII;
	}
	else{
		cptr->Q[ihii] = 1;
	}
}

void Radiation::transferRadiation(const double& dt) const {
	if (nstars > 0) {
		GridCell* srcptr = stars[0].fcausal;
		bool average = true;
		/** Loops over cells in grid causally. */
		for(GridCell* cptr = srcptr; cptr != NULL; cptr = cptr->nextcausal){
		//for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
			double dist2 = 0;
			for (int i = 0; i < gptr->gparams->ND; ++i)
				dist2 += (cptr->xc[i] + stars[0].mod[i] - stars[0].x[i])*(cptr->xc[i] + stars[0].mod[i] - stars[0].x[i]);
			updateColDen(cptr, dist2);
			updateTauSC2(average==false, cptr, dist2);
			if (rparams->SCHEME == IMPLICIT)
				updateTauSC2(average==true, cptr, dist2);
			update_HIIfrac(dt, cptr, 0);
			update_dtau(cptr, 0);
		}
	}
}

void Radiation::transferRadiation(const double& dt, MPIHandler& mpih) const {
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
				part = gptr->boundaries[1];
			}
			int id = 0;
			for (int i = 0; i < (int)part->ghostcells.size(); ++i) {
				for (int j = 0; j < (int)part->ghostcells[i].size(); ++j) {
					for(int ir = 0; ir < NR; ++ir)
						part->ghostcells[i][j]->R[ir] = msgArray[id++];
				}
			}
		}
		transferRadiation(dt);
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

void Radiation::applySrcTerms(const double& dt, const HydroParameters* hparams) {
	if (nstars > 0) {
		for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next) {
			//double HII_old, mu_old, T_old;
			double HII_new, mu_new, E_old, E_new, T_new;
			//HII_old = cptr->U[ihii]/cptr->U[iden];
			HII_new = cptr->Q[ihii];
			//mu_old = 1.0/(HII_old + 1.0);
			mu_new = 1.0/(HII_new + 1.0);
			//T_old = (TMIN + (TMAX-TMIN)*HII_old) / (P_SCALE/RHO_SCALE);
			double ke = 0.0;
			for (int dim = 0; dim < gptr->gparams->ND; ++dim)
			   ke += 0.5*cptr->U[ivel+dim]*cptr->U[ivel+dim]/cptr->U[iden];
			E_old = cptr->U[ipre];
			T_new = (rparams->THI + (rparams->THII-rparams->THI)*HII_new);
			double pre = GAS_CONST*cptr->Q[iden]*scale->RHO*T_new/mu_new;
			E_new = (pre/scale->P)/(hparams->GAMMA-1.0) + ke;

			//////////////
			/** ISOTHERMAL */
			cptr->U[ipre] += E_new - E_old;
			//////////////

			/////////////////////
			/** NOT ISOTHERMAL */
			/*
			double dist2 = 0;
			for (int i = 0; i < gptr->gparams->ND; ++i)
				dist2 += (cptr->xc[i] + stars[0].mod[i] - stars[0].x[i])*(cptr->xc[i] + stars[0].mod[i] - stars[0].x[i])*gptr->dx[i]*gptr->dx[i];
			dist2 *= scale->L*scale->L;
			double f_fuv = Cooling::fluxFUV(0.5*rparams->SOURCE_S*(1.0/scale->T), dist2);
			double av_fuv = Cooling::visualExtinction(5.0e-22/(scale->L*scale->L), cptr->R[icolden]);
			double nh = (cptr->Q[iden] / rparams->H_MASS)*(1.0/(scale->L*scale->L*scale->L));
			*/
			//Collisional ionization cooling.
			//TODO: add effects of collisional ionization cooling.
			/*
			if (rparams->COLL_IONIZATION) {
				temp  = coll_ion_rate(T)*(1.0-P[lv_Hp])*P[lv_Hp]*P[lv_nh];
				R[lv_Hp]   += temp; // coll.ion. to H+ adds to R[i]
				R[lv_eint] -= temp*coll_ion_energy(T); // reduces energy by the amount it took to ionise ion (i-1)
			}
			*/

			//Radiative recombination cooling.
			/*
			cptr->U[ihii] -= cptr->Q[irho]*alphaB(cptr)*cptr->Q[ihii]*cptr->Q[ihii]*(1.0-cptr->Q[ihii])]; // rate [1/s]
			R[lv_Hp] -= temp; // recomb from i to i-1 reduces fraction.
			rate += calcRecombinationEnergy(T)*cptr->Q[ihii]*cptr->Q[ihii]*cptr->Q[iden] / rparams->H_MASS;
			*/

			//Remove energy via cooling.
			//cptr->U[ipre] -= (cptr->R[icool] + Cooling::coolingRate(nh, cptr->Q[ihii], temperature(cptr), f_fuv, av_fuv)/scale->P)*dt;
			/////////////////////

			cptr->U[ihii] = cptr->Q[ihii]*cptr->Q[iden];

		}
	}
}

void Radiation::updateSrcTerms(const double& dt, const HydroParameters* hparams) {
	if (nstars > 0) {
		for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next) {
			//double HII_old, mu_old, T_old;
			double HII_new, mu_new, E_old, E_new, T_new;
			//HII_old = cptr->U[ihii]/cptr->U[iden];
			HII_new = cptr->Q[ihii];
			//mu_old = 1.0/(HII_old + 1.0);
			mu_new = 1.0/(HII_new + 1.0);
			//T_old = (TMIN + (TMAX-TMIN)*HII_old) / (P_SCALE/RHO_SCALE);
			double ke = 0.0;
			for (int dim = 0; dim < gptr->gparams->ND; ++dim)
			   ke += 0.5*cptr->U[ivel+dim]*cptr->U[ivel+dim]/cptr->U[iden];
			E_old = cptr->U[ipre];
			T_new = (rparams->THI + (rparams->THII-rparams->THI)*HII_new);
			double pre = GAS_CONST*cptr->Q[iden]*scale->RHO*T_new/mu_new;
			E_new = (pre/scale->P)/(hparams->GAMMA-1.0) + ke;

			//////////////
			/** ISOTHERMAL */
			cptr->UDOT[ipre] += (E_new-E_old)/dt;
			//////////////

			/////////////////////
			/** NOT ISOTHERMAL */
			/*
			double dist2 = 0;
			for (int i = 0; i < gptr->gparams->ND; ++i)
				dist2 += (cptr->xc[i] + stars[0].mod[i] - stars[0].x[i])*(cptr->xc[i] + stars[0].mod[i] - stars[0].x[i])*gptr->dx[i]*gptr->dx[i];
			dist2 *= scale->L*scale->L;
			double f_fuv = Cooling::fluxFUV(0.5*rparams->SOURCE_S*(1.0/scale->T), dist2);
			double av_fuv = Cooling::visualExtinction(5.0e-22/(scale->L*scale->L), cptr->R[icolden]);
			double nh = (cptr->Q[iden] / rparams->H_MASS)*(1.0/(scale->L*scale->L*scale->L));
			*/
			//Collisional ionization cooling.
			//TODO: add effects of collisional ionization cooling.
			/*
			if (rparams->COLL_IONIZATION) {
				temp  = coll_ion_rate(T)*(1.0-P[lv_Hp])*P[lv_Hp]*P[lv_nh];
				R[lv_Hp]   += temp; // coll.ion. to H+ adds to R[i]
				R[lv_eint] -= temp*coll_ion_energy(T); // reduces energy by the amount it took to ionise ion (i-1)
			}
			*/

			//Radiative recombination cooling.
			/*
			cptr->UDOT[ihii] -= cptr->Q[irho]*alphaB(cptr)*cptr->Q[ihii]*cptr->Q[ihii]*(1.0-cptr->Q[ihii])]; // rate [1/s]
			R[lv_Hp] -= temp; // recomb from i to i-1 reduces fraction.
			rate += calcRecombinationEnergy(T)*cptr->Q[ihii]*cptr->Q[ihii]*cptr->Q[iden] / rparams->H_MASS;
			*/

			//Remove energy via cooling.
			//cptr->UDOT[ipre] -= (cptr->R[icool] + Cooling::coolingRate(nh, cptr->Q[ihii], temperature(cptr), f_fuv, av_fuv)/scale->P)*dt;
			/////////////////////

			cptr->UDOT[ihii] += (cptr->Q[ihii]*cptr->Q[iden] - cptr->U[ihii])/dt;

		}
	}
}

void Radiation::addStar(Star src, MPIHandler& mpih, bool snapToFace[6]) {
	if (src.x[0] >= gptr->fcell->xc[0] && src.x[0] <= gptr->lcell->xc[0]+1) {
		src.fcausal = gptr->locate(src.x[0], src.x[1], src.x[2]);
		src.setCore(mpih.getRank());
	}
	else if (src.x[0] < gptr->fcell->xc[0]) {
		src.fcausal = gptr->locate(gptr->fcell->xc[0], src.x[1], src.x[2]);
		src.setCore(mpih.getRank()-1);
	}
	else if (src.x[0] > gptr->lcell->xc[0]+1) {
		src.fcausal = gptr->locate(gptr->lcell->xc[0], src.x[1], src.x[2]);
		src.setCore(mpih.getRank()+1);
	}
	if (snapToFace[0] && src.x[0] == 0)
		src.mod[0] = 0.5;
	else if (snapToFace[3] && src.x[0] == gptr->gparams->NCELLS[0]-1)
		src.mod[0] = -0.5;
	else if (snapToFace[0] || snapToFace[3]) {
		std::cout << "ERROR: star not in appropriate position for snapping to face - try ";
		std::cout << "placing the star next to the external boundary you would like to snap it to." << std::endl;
		exit(EXIT_FAILURE);
	}
	if (snapToFace[1] && src.x[1] == 0)
		src.mod[1] = 0.5;
	else if (snapToFace[4] && src.x[1] == gptr->gparams->NCELLS[1] -1)
		src.mod[1] = -0.5;
	else if (snapToFace[1] || snapToFace[4]) {
		std::cout << "ERROR: star not in appropriate position for snapping to face - try ";
		std::cout << "placing the star next to the external boundary you would like to snap it to." << std::endl;
		exit(EXIT_FAILURE);
	}
	if (snapToFace[2] && src.x[2] == 0)
			src.mod[2] = 0.5;
	else if (snapToFace[5] && src.x[2] == gptr->gparams->NCELLS[2] -1)
		src.mod[2] = -0.5;
	else if (snapToFace[3] || snapToFace[5]) {
		std::cout << "ERROR: star not in appropriate position for snapping to face - try ";
		std::cout << "placing the star next to the external boundary you would like to snap it to." << std::endl;
		exit(EXIT_FAILURE);
	}
	if (gptr->gparams->ND < 3)
		src.mod[2] = 0;
	if (gptr->gparams->ND < 2)
		src.mod[1] = 0;
	stars.push_back(src);
	nstars++;
	gptr->buildCausal(src.fcausal);
}

bool Radiation::isStar(GridCell* cptr) const {
	for (int i = 0; i < nstars; ++i) {
		if (cptr->xc[0] == stars[i].x[0] && cptr->xc[1] == stars[i].x[1] && cptr->xc[2] == stars[i].x[2]) {
			if (stars[i].mod[0] == 0 && stars[i].mod[1] == 0 && stars[i].mod[2] == 0)
			return true;
		}
	}
	return false;
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
