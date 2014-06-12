/**
 * @file parameters.cpp
 */

#include "parameters.hpp"

#include <stdio.h>

IntegrationParameters::IntegrationParameters(){
	ORDER_S = 1; // Order of reconstruction in cell {CONSTANT, LINEAR, PARABOLIC, etc...}.
	ORDER_T = 1; // Order of accuracy in a time step.
	DT_MAX = 0.001; // Maximum timestep.
}

GridParameters::GridParameters(){
	ND = 1; // No. of spatial dimensions.
	NU = 6; // No. of conserved variables involved in riemann solver.
	NR = 5; // No. of radiation variables not involved in riemann solver.
	NCELLS[0] = 1; // No. of cells along 1st dimension.
	NCELLS[1] = 1; // No. of cells along 2nd dimension.
	NCELLS[2] = 1; // No. of cells along 3rd dimension.
	CORECELLS[0] = 1; // No. of cells along 1st dimension.
	CORECELLS[1] = 1; // No. of cells along 2nd dimension.
	CORECELLS[2] = 1; // No. of cells along 3rd dimension.
	GEOMETRY = CARTESIAN; // Geometry of grid {CARTESIAN, CYLINDRICAL, SPHERICAL}.
	LBCondition[0] = FREE;
	LBCondition[1] = FREE;
	LBCondition[2] = FREE;
	RBCondition[0] = FREE;
	RBCondition[1] = FREE;
	RBCondition[2] = FREE;
}
RadiationParameters::RadiationParameters(){
	K1 = 0.2; // [Mackey 2012 (table 1)] timestep constant.
	K2 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	K3 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	K4 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	P_I_CROSS_SECTION = 6.3E-18; // P.I. cross-sect (cm^2).
	ALPHA_B = 2.59E-13; // Case B radiative recomb rate (cm^3 t^-1 / SCALE).
	TAU_0 = 0.6; // Min. tau in nearest neighbour weights [Mellema et. al. 2006 (A.5)].
	SOURCE_S = 1e51; // Source photon Luminosity (s^-1).
	NHI = 4000; // Initial number density of neutral hydrogen (cm^-3).
	THI = 100; // Temperature fix for fully neutral gas.
	THII = 10000; // Temperature fix for fully ionized gas.
	SCHEME = IMPLICIT; // Ionization fraction integration scheme.
	H_MASS = 1;
}
HydroParameters::HydroParameters(){
	GAMMA = 1.00000001; // Adiabatic gas constant.
	DFLOOR = 0.00001; // Density floor to prevent negative density.
	PFLOOR = 0.00001; // Pressure floor to prevent negative pressure.
}
PrintParameters::PrintParameters(){
	DIR_2D = "tmp/";
	DIR_IF = "tmp/";
	PRINT2D_ON = false;
	PRINTIF_ON = false;
}
Scalings::Scalings(){
	L = 14.6*PC2CM; // cm
	M = 1.0*MO2G; // g
	T = 10000.0*YR2S; // s
	V = L/T;
	RHO = M/(L*L*L);
	P = RHO*V*V;
	E = P*L*L*L;
}
void Scalings::set_LMT(double length, double mass, double time){
	L = length;
	M = mass;
	T = time;
	V = L/T;
	RHO = M/(L*L*L);
	P = RHO*V*V;
	E = P;
}

