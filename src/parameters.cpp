/**
 * @file parameters.cpp
 */

#include "parameters.h"

#include <stdio.h>
#include <iostream>
#include <cmath>

IntegrationParameters::IntegrationParameters(){
	ORDER_S = 1; // Order of reconstruction in cell {CONSTANT, LINEAR, PARABOLIC, etc...}.
	ORDER_T = 1; // Order of accuracy in a time step.
	DT_MAX = 0.001; // Maximum timestep.
}

void IntegrationParameters::print() {
	std::cout << "ORDER_S = " << ORDER_S << "\n";
	std::cout << "ORDER_T = " << ORDER_T << "\n";
	std::cout << "DT_MAX = " << DT_MAX << "\n";
}

GridParameters::GridParameters() {
	ND = 1; // No. of spatial dimensions.
	NCELLS[0] = 1; // No. of cells along 1st dimension.
	NCELLS[1] = 1; // No. of cells along 2nd dimension.
	NCELLS[2] = 1; // No. of cells along 3rd dimension.
	CORECELLS[0] = 1; // No. of cells along 1st dimension.
	CORECELLS[1] = 1; // No. of cells along 2nd dimension.
	CORECELLS[2] = 1; // No. of cells along 3rd dimension.
	GEOMETRY = CARTESIAN; // Geometry of grid {CARTESIAN, CYLINDRICAL, SPHERICAL}.
	SIDE_LENGTH = 1.0;
	LBCondition[0] = FREE;
	LBCondition[1] = FREE;
	LBCondition[2] = FREE;
	RBCondition[0] = FREE;
	RBCondition[1] = FREE;
	RBCondition[2] = FREE;
}

void GridParameters::print() {
	std::cout << "ND = " << ND << "\n";
	std::cout << "NCELLS[0] = " << NCELLS[0] << "\n";
	std::cout << "NCELLS[1] = " << NCELLS[1] << "\n";
	std::cout << "NCELLS[2] = " << NCELLS[2] << "\n";
	std::cout << "CORECELLS[0] = " << CORECELLS[0] << "\n";
	std::cout << "CORECELLS[1] = " << CORECELLS[1] << "\n";
	std::cout << "CORECELLS[2] = " << CORECELLS[2] << "\n";
	std::cout << "GEOMETRY = " << GEOMETRY << "\n";
	std::cout << "SIDE_LENGTH = " << SIDE_LENGTH << "\n";
	std::cout << "LBCondition[0] = " << LBCondition[0] << "\n";
	std::cout << "LBCondition[1] = " << LBCondition[1] << "\n";
	std::cout << "LBCondition[2] = " << LBCondition[2] << "\n";
	std::cout << "RBCondition[0] = " << RBCondition[0] << "\n";
	std::cout << "RBCondition[1] = " << RBCondition[1] << "\n";
	std::cout << "RBCondition[2] = " << RBCondition[2] << "\n";
}

StarParameters::StarParameters() {
	POSITION[0] = 0;
	POSITION[1] = 0;
	POSITION[2] = 0;
	FACE_SNAP[0] = false;
	FACE_SNAP[1] = false;
	FACE_SNAP[2] = false;
	FACE_SNAP[3] = false;
	FACE_SNAP[4] = false;
	FACE_SNAP[5] = false;
	PHOTON_RATE = 0;
	PHOTON_ENERGY = 0;
}

void StarParameters::print() {
	std::cout << "POSITION[0] = " << POSITION[0] << "\n";
	std::cout << "POSITION[1] = " << POSITION[1] << "\n";
	std::cout << "POSITION[2] = " << POSITION[2] << "\n";
	std::cout << "FACE_SNAP[0] = " << FACE_SNAP[0] << "\n";
	std::cout << "FACE_SNAP[1] = " << FACE_SNAP[1] << "\n";
	std::cout << "FACE_SNAP[2] = " << FACE_SNAP[2] << "\n";
	std::cout << "FACE_SNAP[0] = " << FACE_SNAP[3] << "\n";
	std::cout << "FACE_SNAP[1] = " << FACE_SNAP[4] << "\n";
	std::cout << "FACE_SNAP[2] = " << FACE_SNAP[5] << "\n";
	std::cout << "PHOTON_RATE = " << PHOTON_RATE << "\n";
	std::cout << "PHOTON_ENERGY = " << PHOTON_ENERGY << "\n";
}

RadiationParameters::RadiationParameters() {
	K1 = 0.2; // [Mackey 2012 (table 1)] timestep constant.
	K2 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	K3 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	K4 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	P_I_CROSS_SECTION = 6.3E-18; // P.I. cross-sect (cm^2).
	ALPHA_B = 2.59E-13; // Case B radiative recomb rate (cm^3 t^-1 / SCALE).
	TAU_0 = 0.6; // Min. tau in nearest neighbour weights [Mellema et. al. 2006 (A.5)].
	NHI = 4000; // Initial number density of neutral hydrogen (cm^-3).
	THI = 100; // Temperature fix for fully neutral gas.
	THII = 10000; // Temperature fix for fully ionized gas.
	SCHEME = IMPLICIT; // Ionization fraction integration scheme.
	H_MASS = 1;
	SPEC_GAS_CONST = 8.314462e7;
}

void RadiationParameters::print() {
	std::cout << "K1 = " << K1 << "\n";
	std::cout << "K2 = " << K2 << "\n";
	std::cout << "K3 = " << K3 << "\n";
	std::cout << "K4 = " << K4 << "\n";
	std::cout << "P_I_CROSS_SECTION = " << P_I_CROSS_SECTION << "\n";
	std::cout << "ALPHA_B = " << ALPHA_B << "\n";
	std::cout << "TAU_0 = " << TAU_0 << "\n";
	std::cout << "NHI = " << NHI << "\n";
	std::cout << "THI = " << THI << "\n";
	std::cout << "THII = " << THII << "\n";
	std::cout << "SCHEME = " << SCHEME << "\n";
	std::cout << "SPEC_GAS_CONST = " << SPEC_GAS_CONST << "\n";
	for (int i = 0; i < vStarParams.size(); ++i) {
		std::cout << "----STAR[" << i << "]----\n";
		vStarParams[i].print();
	}
}

HydroParameters::HydroParameters(){
	GAMMA = 1.00000001; // Adiabatic gas constant.
	DFLOOR = 0.00001; // Density floor to prevent negative density.
	PFLOOR = 0.00001; // Pressure floor to prevent negative pressure.
}

void HydroParameters::print() {
	std::cout << "GAMMA = " << GAMMA << "\n";
	std::cout << "DFLOOR = " << DFLOOR << "\n";
	std::cout << "PFLOOR = " << PFLOOR << "\n";
}

PrintParameters::PrintParameters(){
	DIR_2D = "tmp/";
	DIR_IF = "tmp/";
	PRINT2D_ON = false;
	PRINTIF_ON = false;
}

void PrintParameters::print() {
	std::cout << "DIR_2D = " << DIR_2D << "\n";
	std::cout << "DIR_IF = " << DIR_IF << "\n";
	std::cout << "PRINT2D_ON = " << PRINT2D_ON << "\n";
	std::cout << "PRINTIF_ON = " << PRINTIF_ON << "\n";
}

Scalings::Scalings(){
	M = 1.0*MO2G; // g
	L = 14.6*PC2CM; // cm
	T = 10000.0*YR2S; // s
	V = L/T;
	RHO = M/(L*L*L);
	P = RHO*V*V;
	E = P*L*L*L;
}
void Scalings::set_mass_length_time(const double& mass, const double& length, const double& time){
	M = mass;
	L = length;
	T = time;
	V = L/T;
	RHO = M/(L*L*L);
	P = RHO*V*V;
	E = P;
}

void Scalings::set_rho_pressure_time(const double& rho, const double& pressure, const double& time){
	T = time;
	M = std::sqrt((pressure/rho)*pressure*pressure)*T*T*T;
	L = std::sqrt(pressure/rho)*T;
	RHO = rho;
	P = pressure;
	E = pressure*L*L*L;
}

double Scalings::toCodeUnits(const double& val, const double& mass_index, const double& length_index, const double& time_index) const {
	return convertCodeUnits(val, mass_index, length_index, time_index, true);
}

double Scalings::fromCodeUnits(const double& val, const double& mass_index, const double& length_index, const double& time_index) const {
	return convertCodeUnits(val, mass_index, length_index, time_index, false);
}

double Scalings::convertCodeUnits(const double& val, const double& mass_index, const double& length_index, const double& time_index,
		const bool& to_code_units) const {
	double MM, LL, TT;
	MM = mass_index >= 0 ? M : 1.0/M;
	MM = to_code_units ? 1.0/MM : MM;
	LL = length_index >= 0 ? L : 1.0/L;
	LL = to_code_units ? 1.0/LL : LL;
	TT = time_index >= 0 ? T : 1.0/T;
	TT = to_code_units ? 1.0/TT : TT;
	int im = (int)std::abs(mass_index);
	int il = (int)std::abs(length_index);
	int it = (int)std::abs(time_index);
	double result = val;
	while (im >= 1 || il >= 1 || it >= 1) {
		if (im-- >= 1) result *= MM;
		if (il-- >= 1) result *= LL;
		if (it-- >= 1) result *= TT;
	}
	double left = std::abs(mass_index) - (int)std::abs(mass_index);
	if (left != 0)
		result *= std::pow(MM, left);
	left = std::abs(length_index) - (int)std::abs(length_index);
	if (left != 0)
		result *= std::pow(LL, left);
	left = std::abs(time_index) - (int)std::abs(time_index);
	if (left != 0)
		result *= std::pow(TT, left);
	return result;
}

void Scalings::print() {
	std::cout << "M = " << M << "\n";
	std::cout << "L = " << L << "\n";
	std::cout << "T = " << T << "\n";
	std::cout << "RHO = " << RHO << "\n";
	std::cout << "P = " << P << "\n";
	std::cout << "E = " << E << "\n";
}
