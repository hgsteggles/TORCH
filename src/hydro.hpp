/* hydro.cpp */

/**
 * Provides the GridCell and GridJoin classes.
 * @file gridcell.hpp
 *
 * @author Harrison Steggles
 * @date 13/01/2014, the first version.
 * @date 16/01/2014, modified code to accomodate new Boundary class.
 */

#ifndef HYDRO_H
#define HYDRO_H

#include <stdlib.h>
#include <math.h>
#include <limits>
#include "boundary.hpp"
#include "grid3d.hpp"
#include "parameters.hpp"
#include "rtmodule.hpp"

/**
 * @class HydroDynamics
 * @brief Contains parameters and methods for integrating the hydrodynamic fluid variables.
 * @version 0.3, 29/01/2014
 */
class HydroDynamics{
public:
	double GAMMA, DFLOOR, PFLOOR, DTMAX;
	HydroDynamics(const HydroParameters& hp);
	void globalWfromU(Grid3D* gptr) const;
	void globalUfromW(Grid3D* gptr) const;
	void globalQfromU(Grid3D* gptr) const;
	void globalUfromQ(Grid3D* gptr) const;
	void UfromQ(double u[], double q[]) const;
	void QfromU(double q[], double u[]) const;
	void FfromU(double f[], double u[], int dim) const;
	double soundSpeed(double pre, double den) const;
	double CFL(Grid3D* gptr) const;
	double fluidStep(Grid3D* gptr) const;
	void calcFluxes(Grid3D* gptr) const;
	void calcBoundaryFluxes(Grid3D* gptr) const;
	void advSolution(double dt, Grid3D* gptr) const;
	void updateBoundaries(Grid3D* gptr) const;
	void HLLC(double U_l[], double F[], double U_r[], int dim) const;
	void reconstruct(Grid3D* gptr) const;
	double av(double a, double b) const;
	void applySrcTerms(double dt, Grid3D* gptr, const Radiation& rad) const;
	void fixSolution(Grid3D* gptr) const;
	void Qisnan(int id, int i, int xc, int yc, int zc, Grid3D* gptr) const;
};

#endif
