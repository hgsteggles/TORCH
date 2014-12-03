/** Provides the Radiation class.
 *
 * @file Radiation.hpp
 *
 * @author Harrison Steggles
 *
 * @date 29/01/2014 - The first version.
 * @date 30/01/2014 - Star passed into addStar can now have any location. Pointer to first GridCell in causal
 * loop (fcausal) is nearest to Star with preference to positive directions.
 * @date 31/01/2014 - Star can only snap to grid vertex at the corner of a GridCell nearest the origin.
 * Coordinates passed in are the grid coordinates of the GridCell. Therefore coordinate transforms should
 * be considered when comparing locations. GridCell::xc on Star x axis is = GridCell::xc + 0.5.
 * @date 31/01/2014 - Fixed bug in first if statement in Radiation::updateTauSC. Now escapes method if Star
 * is located in one of GridCells corners.
 * @date 04/02/2014 - Arguments now passed by const reference when appropriate.
 * @date 13/02/2014 - 28% speedup achieved by adding GridCell member variables that hold cell path length,
 * pointers to nearest neighbours and the weights associated with them instead of calculating them at each
 * time step.
 * @date 11/03/2014 - Bug in transferRadiation(const double, double&, MPIHandler&). Was accessing boundaries[3]
 * instead of boundaries[1].
 * @date 14/04/2014 - Radiation parameters now all held within a RadiationParameters object.
 * @date 14/04/2014 - Function updateColDen added that calculates column density for total hydrogen density
 * (rather than ionized hydrogen as in updateTauSC).
 * @date 14/04/2014 - Function collIonRate added in anticipation of support for collisional ionizations.
 * @date 01/05/2014 - Radiation#getTimeStep is now Radiation#calcTimeStep. More accurately named after what it does.
 * @date 28/05/2014 - New methods added to update source terms instead of directly modifying the conservative fluid variables.
 * @date 13/06/2014 - Renamed PIrate and collIonRate to photoionisationRate and collisionalIonisationRate to make clear
 * what they calculate.
 * @date 02/07/2014 - Radiation#addStar now properly calculates the the core the Star is located in.
 * @date 21/07/2014 - Replaced const reference to int/double with a const copy. The function addStar(Star, bool[]) has
 * been modified so that it only adds the Star to the list of radiating stars. It no longer builds the Star as Stars can
 * now build themselves.
 * @date 24/11/2014 - lots of restructuring and changes to GridCell iteration. Radiation is now an Integrator subclass.
 */

#ifndef RADIATION_H
#define RADIATION_H

#include "Integrator.hpp"
#include "Constants.hpp"
#include "Common.hpp"
#include "SplineData.hpp"

#include <stdlib.h>
#include <vector>
#include <memory>
#include <array>

class GridCell;
class Fluid;
class RadiationParameters;
class HydroParameters;
class Thermodynamics;
class Star;
class StarParameters;
class Converter;

/**
 * @class Radiation
 *
 * @brief Contains parameters and methods for calculating the radiation field around a Star object in a grid.
 *
 * @see RadiationParameters
 * @see Fluid
 * @see Grid
 *
 * @version 0.8, 24/11/2014
 */
class Radiation : public Integrator {
public:
	Radiation();
	~Radiation() {}
	void initialise(std::shared_ptr<Constants> c, RadiationParameters rp);

	void initField(Fluid& fluid) const;

	virtual void preTimeStepCalculations(Fluid& fluid) const;
	virtual double calculateTimeStep(double dt_max, Fluid& fluid) const;
	virtual void integrate(double dt, Fluid& fluid) const;
	virtual void updateSourceTerms(double dt, Fluid& fluid) const;

	double K1 = 0;
	double K2 = 0;
	double K3 = 0;
	double K4 = 0;
	double THI = 0;
	double THII = 0;
	double m_alphaB = 0;
	bool collisions_on = false;
	Coupling coupling = Coupling::OFF;
	double heatingAmplification = 0;
	double massFractionH = 1.0;
	double minX = 0;
	double photoIonCrossSection = 0;
	Scheme scheme = Scheme::IMPLICIT;
	double tau0 = 0;

	std::string printInfo() const;
private:
	std::shared_ptr<Constants> m_consts = nullptr;
	//Hummer (1994) parameters.
	std::unique_ptr<LinearSplineData> m_recombinationHII_CoolingRates;
	std::unique_ptr<LinearSplineData> m_recombinationHII_RecombRates;

	//Initialization methods.
	int getRayPlane(Vec3& xc, Vec3& xs) const;
	// TODO fluid should have this member.
	void calculateNearestNeighbours(Fluid& fluid) const;
	double cellPathLength(Vec3& xc, Vec3& sc, Vec3& dx) const;
	double shellVolume(double ds, double r_sqrd) const;
	void initRecombinationHummer(const Converter& converter);


	//Calculation methods.
	void doric(const double dt, double& HII_avg, double& HII, double Api, double nHII_aB, double nHII_Aci) const;
	double recombinationRateCoefficient(double T) const;
	double recombinationCoolingRate(double nH, double HIIFRAC, double T) const;
	double collisionalIonisationRate(double T) const;
	double photoionisationRate(double nHI, double T, double delT, double shellVol, double photonRate) const;
	double HIIfracRate(double A_pi, double A_ci, double A_rr, double nH, double frac) const;
	double calc_dtau(double nHI, double ds) const;

	//Update methods.
	void updateTauSC(bool average, GridCell& cell, double dist2) const;
	void update_HIIfrac(double dt, GridCell& cell, Fluid& fluid) const;

	//Integration methods.
	void transferRadiation(double dt, Fluid& fluid) const;
	void rayTrace(Fluid& fluid) const;
	void transferRadiation2(double dt, Fluid& fluid) const;

	//Misc. methods.
	bool isStar(const GridCell& cell, const Star& star) const;
};

#endif
