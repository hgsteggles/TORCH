/* io.h */

#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <map>
#include "parameters.hpp"
#include "grid3d.hpp"
#include "rtmodule.hpp"

class InputOutput{
public:
	string DIR_2D, DIR_IF;
	double L_SCALE, M_SCALE, T_SCALE, V_SCALE, RHO_SCALE, P_SCALE, E_SCALE;
	InputOutput(const PrintParameters& rp, const Scalings& sc);
	void print2D(int step, double t, double dt, Grid3D* gptr) const;
	void printIF(GridCell* srcptr, Grid3D* gptr, const Radiation& rad, double t) const;
	void printParams(double runtime, const GridParameters& gpar, const RadiationParameters& rpar, const HydroParameters& hpar, const PrintParameters& ppar, const Scalings& scale) const;
	void fileToMap(const string& myString, std::map<double,double> myMap) const;
};

double getStromRinf(double srcS, double rho, double alpha);
double getStromR(double srcS, double rho, double alpha, double t);
int digits(int x);


#endif
