/* io.cpp */
#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <cmath>
#include "io.hpp"

InputOutput::InputOutput(const PrintParameters& pp, const Scalings& sc) {
	DIR_2D = pp.DIR_2D;
	DIR_IF = pp.DIR_IF;
	PRINT2D_ON = pp.PRINT2D_ON;
	PRINTIF_ON = pp.PRINTIF_ON;
	L_SCALE = sc.L;
	M_SCALE = sc.M;
	T_SCALE = sc.T;
	V_SCALE = sc.V;
	RHO_SCALE = sc.RHO;
	P_SCALE = sc.P;
	E_SCALE = sc.E;
}
void InputOutput::print2D(int step, double t, double dt, Grid3D* gptr, MPIHandler& mpih) const {
	if (PRINT2D_ON) {
		/* creating filename */
		ostringstream os;
		os << DIR_2D << "data2D_";
		os << step << ".txt";
		/* opening new file for appending data */
		ofstream outFile(os.str().c_str(), ios::app);
		if(!outFile)
			cerr << "ERROR: Unable to open " << os.str() << '\n';
		/* writing data to file */
		outFile << setprecision(10) << fixed;
		//outFile << "t" << '\t' << "x" << '\t' << "y" << '\t' << "z" << '\t' << "HIIfrac" << '\t' << "rho" << '\t' << "pre" << '\n';
		if (mpih.getRank() != 0) {
			bool ready;
			mpih.receive(mpih.getRank() - 1, PRINT2D_MSG, ready);
			mpih.wait();
		}
		for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
			outFile << cptr->xc[0];
			if(ND > 1)
				outFile << '\t' << cptr->xc[1];
			if(ND > 2)
				outFile << '\t' << cptr->xc[2];
			outFile << '\t'	<< cptr->Q[iden];
			outFile << '\t' << cptr->Q[ipre];
			outFile << '\t' << cptr->Q[ihii];
			outFile << '\t' << cptr->Q[ivel+0];
			if(ND > 1)
				outFile << '\t' << cptr->Q[ivel+1];
			if(ND > 2)
				outFile << '\t' << cptr->Q[ivel+2];
			outFile << endl;
		}
		if (mpih.getRank() != mpih.nProcessors() - 1) {
			bool ready = true;
			mpih.send(mpih.getRank() + 1, PRINT2D_MSG, ready);
		}
		outFile.close();
	}
}
void InputOutput::printIF(double t, const Radiation& rad, Grid3D* gptr, MPIHandler& mpih) const {
	if (PRINTIF_ON) {
		/* IONIZATION FRONT DETECTION */
		double IF = 0;
		bool found = false;
		if (mpih.getRank() != 0) {
			mpih.receive(mpih.getRank() - 1, PRINTIF_NEXT_MSG, found);
			mpih.wait();
		}
		if(!found) {
			for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->traverse1D(0, 1, cptr)) {
				double frac = cptr->U[ihii]/cptr->U[iden];
				if(cptr->xc[0] >= gptr->fcell->xc[0] && frac < 0.5 && !found) { // IF at HII fraction = 0.5 .
					if(cptr != gptr->fcell){
						if(cptr->left[0] != NULL) {
							double fracLeft = cptr->left[0]->U[ihii]/cptr->left[0]->U[iden];
							double interp = (0.5-fracLeft)/(frac-fracLeft);
							IF = cptr->left[0]->xc[0] - gptr->fcell->xc[0] + interp;
						}
						else
							IF = cptr->xc[0] - gptr->fcell->xc[0];
					}
					found = true;
					break;
				}
			}
			if(mpih.getRank() != mpih.nProcessors()-1)
				mpih.send(mpih.getRank()+1, PRINTIF_NEXT_MSG, found);
			if(mpih.getRank() != 0) {
				mpih.send(0, PRINTIF_FOUND_MSG, found);
				if(found)
					mpih.send(0, PRINTIF_IF_MSG, IF);
			}
			else if(!found) {
				for(int i = 1; i < mpih.nProcessors(); i++) {
					mpih.receive(i, PRINTIF_FOUND_MSG, found);
					mpih.wait();
					if(found) {
						mpih.receive(i, PRINTIF_IF_MSG, IF);
						mpih.wait();
						break;
					}
				}
			}
		}
		if(mpih.getRank() == 0) {
			IF /= gptr->NCELLS[0];
			ofstream outFile("tmp/IF.dat", ios::app);
			double n_H, S, RSinf;
			double cII, RI;
			//double RS, t_rec;
			n_H = rad.NHI / (1.0/(L_SCALE*L_SCALE*L_SCALE));
			//t_rec = (1.0/(n_H*ALPHA));
			double t1 = t;
			S = rad.SOURCE_S;
			RSinf = pow((3.0*S)/(4.0*PI*n_H*n_H*rad.ALPHA_B), 1.0/3.0);
			//RS = pow((1.0-exp(-t1/t_rec)), 1.0/3.0)*RSinf;
			cII = sqrt(GAS_CONST*2.0*rad.TMAX) / V_SCALE;
			RI = RSinf*pow(1+(7*cII*t1)/(4*RSinf), 4.0/7.0);
			double t_s = RSinf/cII;
			double error = 0;
			if (RI != 0.0)
				error = fabs(IF-RI)/RI;
			outFile << t/t_s << '\t';
			outFile << IF*L_SCALE*CM2PC << '\t';
			outFile << RI*L_SCALE*CM2PC << '\t';
			outFile << error << '\n';
			outFile.close();
		}
	}
}
/*
void printParams(double runtime, GridParameters& gpar, RadiationParameters& rpar, HydroParameters& hpar, PrintParameters& ppar, Scalings& scale){
	ofstream outFile("tmp/params.txt");
	outFile << setprecision(8);
	outFile << scientific;
	outFile << setw(15) << left << "dt" << '\t' << p.dt * T_SCALE << '\n';
	outFile << setw(15) << left << "nocells_x" << '\t' << p.nocells[0] << '\n';
	outFile << setw(15) << left << "nocells_y" << '\t' << p.nocells[1] << '\n';
	outFile << setw(15) << left << "nocells_z" << '\t' << p.nocells[2] << '\n';
	outFile << setw(15) << left << "x width" << '\t' << L_SCALE << '\n';
	outFile << setw(15) << left << "geom" << '\t' << p.geom << '\n';
	outFile << setw(15) << left << "scheme" << '\t' << p.scheme << '\n';
	outFile << setw(15) << left << "K1" << '\t' << K1 << '\n';
	outFile << setw(15) << left << "K2" << '\t' << K2 << '\n';
	outFile << setw(15) << left << "K3" << '\t' << K3 << '\n';
	outFile << setw(15) << left << "K4" << '\t' << K4 << '\n';
	outFile << setw(15) << left << "hv" << '\t' << p.hv * E_SCALE << '\n';
	outFile << setw(15) << left << "srcS" << '\t' << SOURCE_S * (1.0/T_SCALE) << '\n';
	outFile << setw(15) << left << "tau_0" << '\t' << TAU_0 << endl;
	outFile << setw(15) << left << "PIcrossSect" << '\t' << PI_CROSS_SECT * (L_SCALE*L_SCALE) << '\n';
	outFile << setw(15) << left << "alphaB" << '\t' << ALPHA * (L_SCALE*L_SCALE*L_SCALE/T_SCALE) << '\n';
	outFile << setw(15) << left << "clumpBL_x" << '\t' << p.clump[0] << '\n';
	outFile << setw(15) << left << "clumpBL_y" << '\t' << p.clump[1] << '\n';
	outFile << setw(15) << left << "clumpW" << '\t' << p.clump[2] << '\n';
	outFile << setw(15) << left << "no_t_steps" << '\t' << p.tstepctr << '\n';
	outFile << setw(15) << left << "runtime" << '\t' << runtime << '\n';
	outFile << setw(15) << left << "units" << '\t' << "cgs" << '\n';
	outFile.close();
}
*/
void InputOutput::fileToMap(const string& myString, std::map<double,double> myMap) const {
	ifstream inFile(myString.c_str(), ios::in);
	if(!inFile)
	  cout << "ERROR: unable to open " << myString << '\n';
	else{
		double a,b;
		while(!inFile.eof() ){
			inFile >> a >> b;
			myMap[a] = b;
		}
	}
  inFile.close();
}
