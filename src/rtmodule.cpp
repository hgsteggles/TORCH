/* rtmodule.c */

#include <math.h>

#include <iostream>
#include <fstream>
#include <cstring>
#include "rtmodule.hpp"

Radiation::Radiation(const RadiationParameters& rp) : nstars(0) {
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

void Radiation::addStar(Star src, Grid3D* gptr, MPIHandler& mpih) {
	stars.push_back(src);
	nstars++;
	if (src.xc[0] >= gptr->fcell->xc[0] && src.xc[0] <= gptr->lcell->xc[0]) {
		stars[nstars-1].fcausal = gptr->locate(src.xc[0], src.xc[1], src.xc[2]);
		stars[nstars-1].setCore(mpih.getRank());
	}
	else if (src.xc[0] < gptr->fcell->xc[0]) {
		stars[nstars-1].fcausal = gptr->locate(gptr->fcell->xc[0], src.xc[1], src.xc[2]);
		stars[nstars-1].setCore(mpih.getRank()-1);
	}
	else if (src.xc[0] > gptr->lcell->xc[0]) {
		stars[nstars-1].fcausal = gptr->locate(gptr->lcell->xc[0], src.xc[1], src.xc[2]);
		stars[nstars-1].setCore(mpih.getRank()+1);
	}
}

int Radiation::getRayPlane(GridCell* cptr, int starID, Grid3D* gptr) const {
	int plane = 0, xc[3], xs[3], dxcs[3], dxcs_max;
	for(int i = 0; i < 3; i++){
		xc[i] = cptr->xc[i];
		xs[i] = stars[starID].xc[i];
		dxcs[i] = fabs(xc[i]-xs[i]);
	}
	dxcs_max = max(dxcs[0], max(dxcs[1], dxcs[2]));
	while(dxcs[plane] != dxcs_max)
		plane++;
	return plane;
}
void Radiation::updateTauSC(bool average, GridCell* cptr, int starID, Grid3D* gptr) const {
	if(cptr->xc[0] == stars[starID].xc[0] && cptr->xc[0] == stars[starID].xc[0] && cptr->xc[0] == stars[starID].xc[0]){
		int plane = getRayPlane(cptr, starID, gptr);
		int icoord = (plane+1)%3, jcoord = (plane+2)%3, kcoord = (plane%3);
		int cpos[3] = {cptr->xc[icoord], cptr->xc[jcoord], cptr->xc[kcoord]};
		int spos[3] = {stars[starID].xc[icoord], stars[starID].xc[jcoord], stars[starID].xc[kcoord]};
		int	dxcs[3] = {cpos[0]-spos[0], cpos[1]-spos[1], cpos[2]-spos[2]};
		int sigma[3] = {abs(dxcs[0]), abs(dxcs[1]), abs(dxcs[2])};
		int LR[3] = {abs(dxcs[0]), abs(dxcs[1]), abs(dxcs[2])};
		
		for(int i = 0; i < 3; i++){
			if(sigma[i] != 0){
				sigma[i] = sigma[i]/dxcs[i];
				LR[i] = LR[i]/dxcs[i];
			}
			else{
				sigma[i] = 1;
				LR[i] = 0;
			}
		}
		double ic[3];
		for(int i = 0; i < 3; i++)
			ic[i] = (double)cpos[i]-0.5*((double)sigma[2]*(double)dxcs[i]/(double)dxcs[2]);
		GridCell* e[4];
		e[0] = gptr->traverse3D(icoord, jcoord, kcoord, 0, 0, -LR[2], cptr);
		e[1] = gptr->traverse3D(icoord, jcoord, kcoord, 0, -LR[1], -LR[2], cptr);
		e[2] = gptr->traverse3D(icoord, jcoord, kcoord, -LR[0], 0, -LR[2], cptr);
		e[3] = gptr->traverse3D(icoord, jcoord, kcoord, -LR[0], -LR[1], -LR[2], cptr);
		if (e[0] == NULL && e[1] == NULL && e[2] == NULL && e[3] == NULL) {
			e[0] = gptr->traverseOverJoins3D(icoord, jcoord, kcoord, 0, 0, -LR[2], cptr);
			e[1] = gptr->traverseOverJoins3D(icoord, jcoord, kcoord, 0, -LR[1], -LR[2], cptr);
			e[2] = gptr->traverseOverJoins3D(icoord, jcoord, kcoord, -LR[0], 0, -LR[2], cptr);
			e[3] = gptr->traverseOverJoins3D(icoord, jcoord, kcoord, -LR[0], -LR[1], -LR[2], cptr);
		}
		double delta[2];
		for(int i = 0; i < 2; i++)
			delta[i] = fabs(2.0*ic[i]-2.0*cpos[i]+sigma[i]);
		double w[4] = {0.0, 0.0, 0.0, 0.0};
		if((dxcs[0] != 0) || (dxcs[1] != 0) || (dxcs[2] != 0))
			w[0] = delta[0]*delta[1];
		if(dxcs[1] != 0)
			w[1] = delta[0]*(1.0-delta[1]);
		if(dxcs[0] != 0)
			w[2] = (1.0-delta[0])*delta[1];
		if((dxcs[0] != 0) && (dxcs[1] != 0))
			w[3] = (1.0-delta[0])*(1.0-delta[1]);
		double tau[4] = {0.0, 0.0, 0.0, 0.0};
		double w_raga[4];
		for(int i = 0; i < 4; i++){
			if(e[i] != NULL){
				if(average == false){
					tau[i] = e[i]->R[itau]+e[i]->R[idtau];
					w_raga[i] = w[i]/max(TAU_0, tau[i]);
				}
				else{
					tau[i] = e[i]->R[itauta]+e[i]->R[idtauta];
					w_raga[i] = w[i]/max(TAU_0, tau[i]);
				}
			}
		}
		double sum_w = w_raga[0]+w_raga[1]+w_raga[2]+w_raga[3];
		double newtau = 0.0;
		for(int i = 0; i < 4; i++){
			if(e[i] != NULL){
				w_raga[i] = w_raga[i]/sum_w;
				newtau += w_raga[i]*tau[i];
			}
		}
		if(average == false)
			cptr->R[itau] = newtau;
		else
			cptr->R[itauta] = newtau;
	}
	else{
		if(average == false)
			cptr->R[itau] = 0;
		else
			cptr->R[itauta] = 0;		
	}
}
double Radiation::temperature(GridCell* cptr) const {
	double mol_mass = 1.0/(cptr->Q[ihii] + 1.0);
	double rho = cptr->Q[iden];
	double pre = cptr->Q[ipre];
	return (pre/rho)*(mol_mass/GAS_CONST);
}
double Radiation::alphaB(GridCell* cptr) const {
	//double T = temperature(cptr);
	//double aB = 0.0;
	double aB = ALPHA_B;//*pow(T/10000, -0.7); // cm^3 s^-1
	return aB;
}
double Radiation::cellPathLength(GridCell* cptr, int starID, Grid3D* gptr) const {
	double ds = 0.0, denom;
	double dx[3] = {gptr->dx[0], gptr->dx[1], gptr->dx[2]};
	int dxcs[3];
	for(int i = 0; i < 3; i++)
		dxcs[i] = cptr->xc[i]-stars[starID].xc[i];
	if(dxcs[2]>=dxcs[1] && dxcs[2]>=dxcs[0])
		denom = dxcs[2];
	else if(dxcs[1]>dxcs[2] && dxcs[1]>=dxcs[0])
		denom = dxcs[1];
	else
		denom = dxcs[0];
	if(denom != 0)
		ds = sqrt(dxcs[0]*dxcs[0]*dx[0]*dx[0]+dxcs[1]*dxcs[1]*dx[1]*dx[1]+dxcs[2]*dxcs[2]*dx[2]*dx[2])/denom;
	return ds;
}
double Radiation::shellVolume(GridCell* cptr, int starID, Grid3D* gptr) const {
	double r_sqrd = 0;
	double ds = cellPathLength(cptr, starID, gptr);
	double dx[3] = {gptr->dx[0], gptr->dx[1], gptr->dx[2]};
	for(int i=0; i<3; i++)
		r_sqrd += (cptr->xc[i]-stars[starID].xc[i])*(cptr->xc[i]-stars[starID].xc[i])*dx[i]*dx[i];
	double r_min = sqrt(r_sqrd) - 0.5*ds;
	double r_max = sqrt(r_sqrd) + 0.5*ds;
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
		if(r_sqrd <= 0.1)
			return 4.0*PI*(r_max*r_max*r_max - r_min*r_min*r_min)/3.0;
		else
			return 4.0*PI*r_sqrd*ds;
}
double Radiation::PIrate(double frac, double T, double delT, GridCell* cptr, int starID, Grid3D* gptr) const {
	double vol = shellVolume(cptr, starID, gptr);
	double n_HI = (1.0-frac)*cptr->Q[iden] / H_MASS;
	double A_pi;
	double S = SOURCE_S;// / (1.0/T_SCALE);
	if(frac >= 1 || n_HI == 0 || vol == 0)
		A_pi = 0.0;
	else
		A_pi = S*exp(-T)*(1.0-exp(-delT))/(n_HI*vol);
	return A_pi;
}
double Radiation::HIIfracDot(double A_pi, double frac, GridCell* cptr) const {
	//if(COLLISIONS)
		//A_ci = p.CIlookup[T];
	double nH = cptr->Q[iden] / H_MASS;
	return (1.0-frac)*A_pi - (frac*frac*alphaB(cptr)*nH); // + (1-x)*(x*n_H*A_ci);
}
void Radiation::update_dtau(GridCell* cptr, int starID, Grid3D* gptr) const {
	if(!isStar(cptr)){
		double ds = cellPathLength(cptr, starID, gptr);
		double n_H = cptr->Q[iden] / H_MASS;
		cptr->R[idtau] = (1.0-cptr->Q[ihii])*n_H*P_I_CROSS_SECTION*ds;
		cptr->R[idtauta] = (1.0-cptr->R[ihiita])*n_H*P_I_CROSS_SECTION*ds;
		bool test = false;
		if (cptr->xc[0] == 2 && cptr->xc[1] == 1 && test) {
			cerr << "ds = " << ds << endl;
			cerr << "n_H = " << n_H << endl;
			cerr << "piCrSc = " << P_I_CROSS_SECTION << endl;
			cerr << "HII = " << cptr->Q[ihii] << endl;
			cerr << "dtau = " << cptr->R[idtau] << endl;
		}
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
void Radiation::doric(double dt, double& HII, double& HII_avg, double A_pi, GridCell* cptr) const {
	double HII_old, n_HII_avg, xeq, inv_ti;
	double aB = alphaB(cptr);
	double epsilon = 1.0e-14;
	//double frac_avg_old;
	HII_old = HII;
	//frac_avg_old = frac_avg;
	n_HII_avg = HII_avg*cptr->Q[iden] / H_MASS;
	if(n_HII_avg == 0.0)
		xeq = 1.0;
	else
		xeq = A_pi/(A_pi + n_HII_avg*aB);
	inv_ti = A_pi + n_HII_avg*aB;
	HII = xeq + (HII_old-xeq)*exp(-dt*inv_ti);
	HII = max(min(HII, 1.0), 0.0);
	if(1.0 - HII < epsilon)
		HII = 1.0 - epsilon;
	if(dt*inv_ti < 1.0e-8)
		HII_avg = HII_old;
	else{
		HII_avg = xeq + (HII_old-xeq)*(1.0-exp(-dt*inv_ti))/(dt*inv_ti);
		HII_avg = max(min(HII_avg, 1.0), 0.0);
	}
	if(1.0 - HII_avg < epsilon)
		HII_avg = 1.0 - epsilon;
}
void Radiation::update_HIIfrac(double dt, GridCell* cptr, int starID, Grid3D* gptr) const {
	int scheme = 0;
	if(!isStar(cptr)){
		if(scheme == 0){
			double convergence2 = 1.0e-3;
			double convergence_frac = 1.0e-5;
			double HII;
			double HII_avg_old, HII_avg, dtau_avg, tau_avg, A_pi, ds, n_H;
			int niter = 0;
			bool converged = false;
			ds = cellPathLength(cptr, starID, gptr);
			n_H = cptr->Q[iden] / H_MASS;
			HII = cptr->Q[ihii];
			//HII_avg = cptr->R[ihiita];
			tau_avg = cptr->R[itauta];
			HII_avg = HII;
			//tau_avg = cptr->R[itau];
			while(!converged){
				niter++;
				HII_avg_old = HII_avg;
				HII = cptr->Q[ihii];
				dtau_avg = (1.0-HII_avg)*n_H*P_I_CROSS_SECTION*ds;
				A_pi = PIrate(HII_avg, tau_avg, dtau_avg, cptr, starID, gptr);
				doric(dt, HII, HII_avg, A_pi, cptr);
				if((fabs((HII_avg-HII_avg_old)/HII_avg) < convergence2 || (HII_avg < convergence_frac)))
					converged = true;
				if(niter > 100000){
					std::cout << "ERROR: implicit method not converging." << endl;
					cerr << "cptr->Q[iden] = " << cptr->Q[ihii] << endl;
					cerr << "cptr->Q[ipre] = " << cptr->Q[ihii] << endl;
					cerr << "cptr->Q[ivel+0] = " << cptr->Q[ihii] << endl;
					cerr << "cptr->Q[ihii] = " << cptr->Q[ihii] << endl;
					cerr << "HII_avg_old = " << HII_avg_old << endl;
					cerr << "dtau_avg = " << dtau_avg << endl;
					cerr << "A_pi = " << A_pi << endl;
					cerr << "HII = " << HII << endl;
					cerr << "HII_avg = " << HII_avg_old << endl;
					cerr << "xc[0] = " << cptr->xc[0] << endl;
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
			A_pi = PIrate(HII, tau, dtau, cptr, starID, gptr);
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
double Radiation::getTimeStep(double dt_dyn, Grid3D* gptr) const {
	double dt = dt_dyn, dt1, dt2, dt3, dt4;
	if (nstars > 0) {
		GridCell* srcptr = stars[0].fcausal;
		for(GridCell* cptr = srcptr; cptr != NULL; cptr = gptr->nextCausal(cptr, srcptr)){
			dt1 = dt_dyn;
			dt2 = dt_dyn;
			dt3 = dt_dyn;
			dt4 = dt_dyn;
			double A_pi = PIrate(cptr->Q[ihii], cptr->R[itau], cptr->R[idtau], cptr, 0, gptr);
			double fracRate = HIIfracDot(A_pi, cptr->Q[ihii], cptr);
			if(K1 != 0.0)
				dt1 = K1*1.0/(cptr->Q[iden]*alphaB(cptr)/H_MASS);
			if(K2 != 0.0)
				dt2 = dt_dyn;
			if(K3 != 0.0 && fracRate != 0.0)
				dt3 = K3*max(0.05, 1.0 - cptr->Q[ihii])/fabs(fracRate);
			if(K4 != 0.0 && fracRate != 0.0){
				dt4 = K4*1.0/fabs(fracRate); // timestep criterion [Mackey 2012].
				if(fracRate == 0.0)
					cerr << "fracRate = 0" << endl;
			}
			dt = min(dt, min(dt_dyn, min(dt1, min(dt2, min(dt3, dt4)))));
			if(dt == 0.0){
				cerr << dt << '\t' << dt1 << '\t' << dt_dyn << '\t' << cptr->xc[1] << '\n';
				break;
			}
		}
	}
	return dt;
}
void Radiation::transferRadiation(double dt, double& IF, Grid3D* gptr) const {
	if (nstars > 0) {
		GridCell* srcptr = stars[0].fcausal;
		bool average = true;
		/** LOOPS OVER CELLS IN GRID CAUSALLY.
		 *
		 */
		for(GridCell* cptr = srcptr; cptr != NULL; cptr = gptr->nextCausal(cptr, srcptr)){
			updateTauSC(average==false, cptr, 0, gptr);
			updateTauSC(average==true, cptr, 0, gptr);
			update_HIIfrac(dt, cptr, 0, gptr);
			update_dtau(cptr, 0, gptr);
		}
	}
}

void Radiation::transferRadiation(double dt, double& IF, Grid3D* gptr, MPIHandler& mpih) const {
	if (nstars > 0) {
		vector<double> msgVec(gptr->boundaries[0]->ghostcells.size()*gptr->boundaries[0]->ghostcells[0].size()*NR, 0);
		if(stars[0].core != mpih.getRank()) {
			Boundary* part = NULL;
			if(stars[0].core < mpih.getRank()) {
				mpih.receive(mpih.getRank()-1, PARTITION_MSG, msgVec);
				mpih.wait();
				part = gptr->boundaries[0];
			}
			else {
				mpih.receive(mpih.getRank()+1, PARTITION_MSG, msgVec);
				mpih.wait();
				part = gptr->boundaries[3];
			}
			int id = 0;
			for (int i = 0; i < (int)part->ghostcells.size(); i++) {
				for (int j = 0; j < (int)part->ghostcells[i].size(); j++) {
					for(int ir = 0; ir < NR; ir++) {
						part->ghostcells[i][j]->R[ir] = msgVec[id];
						id++;
					}
				}
			}
		}
		transferRadiation(dt, IF, gptr);
		int face1 = 0, face2 = 3;
		if(mpih.getRank() == 0 || stars[0].core > mpih.getRank())
			face1 = -1;
		if(mpih.getRank() == mpih.nProcessors()-1 || stars[0].core < mpih.getRank())
			face2 = -1;
		for(int i = 0; i < gptr->boundaries.size(); i++) {
			if(gptr->boundaries[i]->face == face1 || gptr->boundaries[i]->face == face2) {
				Partition* part = (Partition*)gptr->boundaries[i];
				int dim = part->face%3;
				int id = 0;
				for (int i = 0; i < (int)part->ghostcells.size(); i++) {
					for (int j = 0; j < (int)part->ghostcells[i].size(); j++) {
						GridCell* cptr = NULL;
						if (part->face < 3)
							cptr = part->ghostcells[i][j]->rjoin[dim]->rcell;
						else
							cptr = part->ghostcells[i][j]->ljoin[dim]->lcell;
						for(int ir = 0; ir < NR; ir++) {
							msgVec[id] = cptr->R[ir];
							id++;
						}
					}
				}
				mpih.send(part->destination, PARTITION_MSG, msgVec);
			}
		}
	}
}

void Radiation::rayTrace(Grid3D* gptr) const {
	static bool time = 0;
	if(time == 0 && nstars > 0){
		GridCell* srcptr = stars[0].fcausal;
		for(GridCell* cptr = srcptr; cptr != NULL; cptr = gptr->nextCausal(cptr, srcptr)){
			if(cptr != srcptr)
				update_dtau(cptr, 0, gptr);
		}
		time++;
	}
}

bool Radiation::isStar(GridCell* cptr) const {
	for (int i = 0; i < nstars; i++) {
		if (cptr->xc[0] == stars[i].xc[0]) {
			if (cptr->xc[1] == stars[i].xc[1]) {
				if (cptr->xc[2] == stars[i].xc[2])
					return true;
			}
		}
	}
	return false;
}
