/**
 * Provides the various Parameter classes.
 * @file parameters.h
 *
 * @author Harrison Steggles
 *
 * @date 13/01/2014 - the first version.
 * @date 31/01/2014 - removed ND and replaced all occurrences of it within CFD files with Grid3D::ND.
 * @date 01/05/2014 - IntegrationParameters class added.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "constants.h"

#include <string>
#include <vector>

class IntegrationParameters{
public:
	IntegrationParameters();
	int ORDER_S, ORDER_T;
	double DT_MAX;
	void print();
};

class GridParameters {
public:
	GridParameters();
	int ND; //!< Number of dimensions.
	int NCELLS[3]; //!< Array holding the number of grid cells along each dimension.
	int CORECELLS[3]; //!< Array holding the number of grid cells along in each dimension in each processor.
	double SIDE_LENGTH; //!< The side length of the simulation line/square/cube.
	Condition LBCondition[3]; //!< Array of left boundary conditions for each dimension.
	Condition RBCondition[3]; //!< Array of right boundary conditions for each dimension.
	Geometry GEOMETRY; //!< The geometry of the grid [CARTESIAN, CYLINDRICAL, SPHERICAL].
	void print();
};

class StarParameters {
public:
	StarParameters();
	int POSITION[3];
	bool FACE_SNAP[6];
	double PHOTON_ENERGY;
	double PHOTON_RATE;
	void print();
};

class RadiationParameters {
public:
	RadiationParameters();
	double K1, K2, K3, K4;
	double P_I_CROSS_SECTION;
	double ALPHA_B;
	double TAU_0;
	double NHI;
	double THI;
	double THII;
	double SCHEME;
	double H_MASS;
	double SPEC_GAS_CONST;
	std::vector<StarParameters> vStarParams;
	void print();
};

class HydroParameters{
public:
	HydroParameters();
	double GAMMA, DFLOOR, PFLOOR;
	void print();
};
class PrintParameters{
public:
	PrintParameters();
	std::string DIR_2D, DIR_IF;
	bool PRINT2D_ON, PRINTIF_ON;
	void print();
};
class Scalings{
public:
	Scalings();
	void set_mass_length_time(const double& mass, const double& length, const double& time);
	void set_rho_pressure_time(const double& rho, const double& pressure, const double& time);
	double toCodeUnits(const double& val, const double& mass_index, const double& length_index, const double& time_index) const;
	double fromCodeUnits(const double& val, const double& mass_index, const double& length_index, const double& time_index) const;
	void print();
private:
	double M, L, T, V, RHO, P, E;
	double convertCodeUnits(const double& val, const double& mass_index, const double& length_index, const double& time_index, const bool& from) const;
};

#endif
