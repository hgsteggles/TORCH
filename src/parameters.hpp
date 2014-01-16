/* parameters.h */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdio.h>
#include <string>
#include <map>
using namespace std;

enum Geometry {CARTESIAN, CYLINDRICAL, SPHERICAL};
enum Condition {FREE, REFLECTING, OUTFLOW, INFLOW};
enum Scheme {IMPLICIT, EXPLICIT};

class GridParameters{
public:
	GridParameters();
	int ND, NU, NR, NCELLS[3], ORDER_S, ORDER_T;
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
	string DIR_2D, DIR_IF;
};
class Scalings{
public:
	Scalings();
	void set_LMT(double L, double M, double T);
	double L, M, T, V, RHO, P, E;
};

const int ND = 1;
const int NU = 6;
const int iden = 0;
const int ivel = 1;
const int ipre = 4;
const int ihii = 5;
const int NR = 5;
const int ihiita = 0;
const int itau = 1;
const int itauta = 2;
const int idtau = 3;
const int idtauta = 4;

extern bool SWEEPX;
/* MICROPHYSICS */
const bool RTCOUPLING = true;
const bool MONOCHROME = true;
const bool RECOMBINATIONS = true;
const bool COLLISIONS = false;
const bool COOLING = false;
/* CONVERSIONS */
const double YR2S = (3.15569e7);
const double S2YR = (1.0/YR2S);
const double PC2CM = (3.09e18);
const double CM2PC = (1.0/PC2CM);
const double MO2G = (2e33);
const double EV2ERG = (1.60217646e-12);
/* CONSTANTS */
const double PI = 3.14159265359;
const double GAS_CONST = 8.314462E7; // ergs mol^-1 K*-1
const double H_MASS_G = (1.674e-24); // g

#endif
