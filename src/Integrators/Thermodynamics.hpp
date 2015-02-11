/** Provides the Thermodynamics class.
 *
 * @file Thermodynamics.hpp
 *
 * @author Harrison Steggles
 *
 * @date 12/06/2014 - The first version.
 * @date 13/06/2014 - Fixed recombinationHII() and added Boltzmann's constant to Thermodynamics as a member variable.
 * @date 02/07/2014 - Everything is in code units now.
 * @date 02/07/2014 - Added method for filling heating arrays for printing.
 * @date 02/07/2014 - Fixed heat flux calculation. Rays are traced across the grid to update column densities.
 * @date 24/11/2014 - Lots of restructuring and changes to GridCell iteration. Thermodynamics is now an Integrator subclass.
 * @date 28/11/2014 - Fixed raytrace bug.
 */

#ifndef THERMODYNAMICS_HPP_
#define THERMODYNAMICS_HPP_

#include "Integrator.hpp"
#include "SplineData.hpp"

#include <vector>
#include <cmath>
#include <memory>

class Converter;
class Star;
class Fluid;
class GridCell;
class MPIW;
class ThermoParameters;
class Constants;

/**
 * @class Thermodynamics
 *
 * @brief Contains parameters and methods for calculating heating and cooling due to atomic processes in the gas around an ionising star.
 *
 * @version 0.8, 24/11/2014
 */
class Thermodynamics : public Integrator {
public:
	Thermodynamics();
	~Thermodynamics() {};

	void initialise(std::shared_ptr<Constants> c, ThermoParameters tp);

	virtual void preTimeStepCalculations(Fluid& fluid) const;
	virtual double calculateTimeStep(double dt_max, Fluid& fluid) const;
	virtual void integrate(double dt, Fluid& fluid) const;
	virtual void updateSourceTerms(double dt, Fluid& fluid) const;

	void fillHeatingArrays(const Fluid& fluid);
private:
	void initCollisionalExcitationHI(const Converter& scale);
	void initRecombinationHII(const Converter& scale);
	double fluxFUV(const double Q_FUV, const double dist_sqrd) const;
	double collisionalExcitationHI(const double nH, const double HIIFRAC, const double T) const;
	double recombinationHII(const double nH, const double HIIFRAC, const double T) const;
	double ionisedMetalLineCooling(const double ne, const double T) const;
	double neutralMetalLineCooling(const double ne, const double nn, const double T) const;
	double collisionalIonisationEquilibriumCooling(const double ne, const double T) const;
	double neutralMolecularLineCooling(const double nH, const double HIIFRAC, const double T) const;
	double farUltraVioletHeating(const double nH, const double Av_FUV, const double F_FUV) const;
	double infraRedHeating(const double nH, const double Av_FUV, const double F_FUV) const;
	double cosmicRayHeating(const double nH) const;
	double softLanding(const double rate, const double T) const;

	void updateColDen(GridCell& cell, const double dist2) const;
	void rayTrace(const Fluid& fluid) const;

	std::shared_ptr<Constants> m_consts = nullptr;

	bool m_isSubcycling = false;
	double m_thermoHII_Switch = 0;
	double m_heatingAmplification = 1.0; //!< Heating amplification/reduction hack.
	double m_massFractionH = 0;

	double m_z0 = 0; //!< ISM oxygen abundance.
	double m_excessEnergy = 0; //!< Excess photon energy [erg].
	double m_T1 = 0;
	double m_T2 = 0;
	double m_T3 = 0;
	double m_T4 = 0;
	double m_imlc = 0;
	double m_nmlc = 0;
	double m_cxhi_damp = 0;
	double m_ciec = 0;
	double m_ciec_minT = 0;
	double m_n0 = 0;
	double m_nmc = 0;
	double m_fuvh_a = 0;
	double m_fuvh_b = 0;
	double m_fuvh_c = 0;
	double m_irh_a = 0;
	double m_irh_b = 0;
	double m_crh = 0;
	double m_T_min = 0;
	double m_T_soft = 0;

	std::vector<double> cxhi_T;
	std::vector<double> cxhi_rate;
	std::vector<double> cxhi_rate2;
	double cxhi_minSlope = 0, cxhi_maxSlope = 0;

	std::vector<double> hr_T;
	std::vector<double> hr_rate;
	std::vector<double> hr_rate2;
	double hr_minSlope = 0, hr_maxSlope = 0;

	std::unique_ptr<LogSplineData> m_collisionalExcitationHI_CoolingRates;
	std::unique_ptr<LinearSplineData> m_recombinationHII_CoolingRates;
};

#endif // THERMODYNAMICS_HPP_
