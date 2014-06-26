/* hydro.cpp */

/**
 * Provides the GridCell and GridJoin classes.
 * @file gridcell.h
 *
 * @author Harrison Steggles
 * @date 13/01/2014 - the first version.
 * @date 16/01/2014 - modified code to accomodate new Boundary class.
 * @date 04/01/2014 - arguments now passed by const reference when appropriate.
 * @date 04/01/2014 - removed unnecessary variable declarations.
 * @date 21/02/2014 - globalUfromW and globalWfromU now uses memcpy rather than copying each element in a for loop.
 * @date 14/04/2014 - Hydro parameters now all held within a HydroParameters object.
 * @date 14/04/2014 - the function HLL has been added to solve the Riemann problem using the HLL scheme.
 * @date 14/04/2014 - new Riemann solver implemented. The function RotatedHLLC is now used to solve the Riemann
 * problem.
 * @date 01/05/2014 - calcFluxes now takes and integer argument: order, rather than using the fixed ORDER_S parameter in Grid3D.
 * @date 01/05/2014 - CFL now takes dt_max as argument instead of setting it from HydroParameters. It doesn't belong in HydroParameters.
 * @date 28/05/2014 - new methods added to update source terms instead of directly modifying the conservative fluid variables.
 */

#ifndef HYDRO_H
#define HYDRO_H

class Grid3D;
class HydroParameters;
class Radiation;
/**
 * @class HydroDynamics
 *
 * @brief Contains parameters and methods for integrating the hydrodynamic fluid variables.
 *
 * @version 0.7, 13/06/2014
 */
class HydroDynamics{
public:
	HydroParameters& hparams;
	Grid3D& grid;
	HydroDynamics(HydroParameters& hp, Grid3D& g3d);
	//Conversion methods.
	void globalWfromU() const;
	void globalUfromW() const;
	void globalQfromU() const;
	void globalUfromQ() const;
	void UfromQ(double u[], double q[]) const;
	void QfromU(double q[], double u[]) const;

	//Calculation methods.
	void reconstruct() const;
	void calcFluxes(const int& order) const;
	double CFL(const double& dt_max) const;

	//Integration methods.
	void updateSrcTerms() const;
	void updateBoundaries() const;
	void advSolution(const double& dt);
	void applySrcTerms(const double& dt) const;
	void fixSolution() const;

private:
	//Calculation methods.
	void FfromU(double f[], double u[], const int& dim) const;
	double soundSpeed(const double& pre, const double& den) const;
	void HLLC(double U_l[], double F[], double U_r[], const int& dim) const;
	void HLL(double Q_l[], double F[], double Q_r[], const int& dim) const;
	void RotatedHLLC(double Q_l[], double F[], double Q_r[], const int& dim) const;

	//Misc. methods.
	double av(const double& a, const double& b) const;
	void Qisnan(const int& id, const int& i, const int& xc, const int& yc, const int& zc) const;
};

#endif
