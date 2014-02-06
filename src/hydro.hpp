/* hydro.cpp */

/**
 * Provides the GridCell and GridJoin classes.
 * @file gridcell.hpp
 *
 * @author Harrison Steggles
 * @date 13/01/2014, the first version.
 * @date 16/01/2014, modified code to accomodate new Boundary class.
 * @date 04/01/2014, arguments now passed by const reference when appropriate.
 * @date 04/01/2014, removed unnecessary variable declarations.
 */

#ifndef HYDRO_H
#define HYDRO_H

class Grid3D;
class HydroParameters;
class Radiation;
/**
 * @class HydroDynamics
 * @brief Contains parameters and methods for integrating the hydrodynamic fluid variables.
 * @version 0.3, 29/01/2014
 */
class HydroDynamics{
public:
	Grid3D* gptr;
	double GAMMA, DFLOOR, PFLOOR, DTMAX;
	HydroDynamics(const HydroParameters& hp, Grid3D* grid);
	void globalWfromU() const;
	void globalUfromW() const;
	void globalQfromU() const;
	void globalUfromQ() const;
	void UfromQ(double u[], double q[]) const;
	void QfromU(double q[], double u[]) const;
	void FfromU(double f[], double u[], const int& dim) const;
	double soundSpeed(const double& pre, const double& den) const;
	double CFL() const;
	double fluidStep() const;
	void calcFluxes() const;
	void calcBoundaryFluxes() const;
	void advSolution(const double& dt) const;
	void updateBoundaries() const;
	void HLLC(double U_l[], double F[], double U_r[], const int& dim) const;
	void reconstruct() const;
	double av(const double& a, const double& b) const;
	void applySrcTerms(const double& dt, const Radiation* rad) const;
	void fixSolution() const;
	void Qisnan(const int& id, const int& i, const int& xc, const int& yc, const int& zc) const;
};

#endif
