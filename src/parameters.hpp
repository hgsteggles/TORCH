/**
 * Provides the various Parameter classes.
 * @file parameters.hpp
 *
 * @author Harrison Steggles
 *
 * @date 13/01/2014, the first version.
 * @date 31/01/2014, removed ND and replaced all occurrences of it within CFD files with Grid3D::ND.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "constants.hpp"

#include <string>

class GridParameters{
public:
	GridParameters();
	int ND, NU, NR, NCELLS[3], ORDER_S, ORDER_T;
	Condition LBCondition[3];
	Condition RBCondition[3];
	Geometry GEOMETRY;
};
class RadiationParameters{
public:
	RadiationParameters();
	double K1, K2, K3, K4, P_I_CROSS_SECTION, ALPHA_B, TAU_0, SOURCE_S, NHI, TMIN, TMAX, SCHEME, H_MASS;
};
class HydroParameters{
public:
	HydroParameters();
	double GAMMA, DFLOOR, PFLOOR, DTMAX;
};
class PrintParameters{
public:
	PrintParameters();
	std::string DIR_2D, DIR_IF;
	bool PRINT2D_ON, PRINTIF_ON;
};
class Scalings{
public:
	Scalings();
	void set_LMT(double L, double M, double T);
	double L, M, T, V, RHO, P, E;
};

#endif
