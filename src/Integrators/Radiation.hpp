/** Provides the Radiation class.
 *
 * @file Radiation.hpp
 *
 * @author Harrison Steggles
 */

#ifndef RADIATION_HPP_
#define RADIATION_HPP_

#include <memory>
#include <string>

#include "Torch/Common.hpp"
#include "Torch/Constants.hpp"
#include "Integrator.hpp"
#include "SplineData.hpp"

class GridCell;
class Fluid;
class RadiationParameters;
class Star;
class StarParameters;
class Converter;

/**
 * @class Radiation
 *
 * @brief Contains parameters and methods for calculating the radiation field around a Star object in a grid. Uses Mellema (2005) C^2-Ray scheme.
 *
 * @see Fluid
 * @see Grid
 * @see Star
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
	std::unique_ptr<LinearSplineData> m_recombinationHII_CoolingRates = nullptr; //!< Hummer (1994) hydrogen recombination cooling rates.
	std::unique_ptr<LinearSplineData> m_recombinationHII_RecombRates = nullptr; //!< Hummer (1994) hydrogen recombination rates.

	// Initialisation methods.
	int getRayPlane(Vec3& xc, Vec3& xs) const;
	double cellPathLength(Vec3& xc, Vec3& sc, Vec3& dx) const;
	double shellVolume(double ds, double r_sqrd) const;
	void initRecombinationHummer(const Converter& converter);


	// Calculation methods.
	void doric(const double dt, double& HII_avg, double& HII, double Api, double nHII_aB, double nHII_Aci) const;
	double recombinationRateCoefficient(double T) const;
	double recombinationCoolingRate(double nH, double HIIFRAC, double T) const;
	double collisionalIonisationRate(double T) const;
	double photoionisationRate(double nHI, double T, double delT, double shellVol, double photonRate) const;
	double HIIfracRate(double A_pi, double A_ci, double A_rr, double nH, double frac) const;
	double calc_dtau(double nHI, double ds) const;

	// Update methods.
	void updateTauSC(bool average, GridCell& cell, Fluid& fluid, double dist2) const;
	void update_HIIfrac(double dt, GridCell& cell, Fluid& fluid) const;

	// Integration methods.
	void transferRadiation(double dt, Fluid& fluid) const;
	void rayTrace(Fluid& fluid) const;
	void transferRadiation2(double dt, Fluid& fluid) const;

	// Misc. methods.
	bool isStar(const GridCell& cell, const Star& star) const;
};

#endif // RADIATION_HPP_
