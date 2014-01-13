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
	L_SCALE = sc.L;
	M_SCALE = sc.M;
	T_SCALE = sc.T;
	V_SCALE = sc.V;
	RHO_SCALE = sc.RHO;
	P_SCALE = sc.P;
	E_SCALE = sc.E;
}
////// digits(x) returns number of digits of an int (for naming
////// files)
int digits(int x){
	int N = 0;
	double y = (double)x;
	while(y >= 1.0){
		y /= 10.0;
		N++;
	}
	if(x == 0)
		N = 1;
	return N;
}
////// print2DToFile prints the cell coordinates in the 1st two
////// columns and the HII fraction of that cell in the 3rd column.
void InputOutput::print2D(int step, double t, double dt, Grid3D* gptr) const {
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
	outFile.close();
}
void InputOutput::printIF(GridCell* srcptr, Grid3D* gptr, const Radiation& rad, double t) const {
	/* IONIZATION FRONT DETECTION */
	double IF = 0;
	bool found = false;
	for(GridCell* cptr = srcptr; cptr != NULL; cptr = gptr->nextCausal(cptr, srcptr)){
		double frac = cptr->U[ihii]/cptr->U[iden];
		if(cptr->xc[0] >= srcptr->xc[0] && frac < 0.5 && !found) { // IF at HII fraction = 0.5 .
			if(cptr != srcptr){
				double fracLeft = cptr->left[0]->U[ihii]/cptr->left[0]->U[iden];
				double interp = (0.5-fracLeft)/(frac-fracLeft);
				IF = cptr->left[0]->xc[0] - srcptr->xc[0] + interp;
			}
			found = true;
			break;
		}
	}
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

////// fileToMap looks up a file and creates a look up table
////// from it.
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
