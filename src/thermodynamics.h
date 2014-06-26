/**
 * Provides the Thermodynamics class.
 *
 * @file thermodynamics.h
 *
 * @author Harrison Steggles
 *
 * @date 12/06/2014 - The first version.
 * @date 13/06/2014 - Fixed recombinationHII() and added Boltzmann's constant to Thermodynamics as a member variable.
 */

#include <vector>
#include <cmath>

class Scalings;
class Radiation;
class Grid3D;

/**
 * @class Thermodynamics
 *
 * @brief Contains parameters and methods for calculating heating and cooling due to atomic processes in the gas around an ionizing star.
 *
 * @version 0.7, 12/06/2014
 */
class Thermodynamics {
public:
	Thermodynamics(const Scalings& scale, Grid3D& g3d);
	double fluxFUV(const double& Q_FUV, const double dist_sqrd);
	double calcTimeStep(const double& dt_dyn);
	void calcHeatFlux(const Radiation& rad);
	double updateSrcTerms();
private:
	Grid3D& grid;
	double spec_gas_const; //!< specific gas constant [g-1.erg.mol-1.K-1]
	double kB; //!< Boltzmann's constant [erg.K-1].
	double sigmaV; //!< Dust extinction cross-section [cm2 H-1].
	double z0; //!< ISM oxygen abundance.
	double excessEnergy; //!< Excess photon energy [erg].
	double T1;
	double T2;
	double T3;
	double T4;
	double imlc;
	double nmlc;
	double cxhi_damp;
	double ciec;
	double ciec_minT;
	double n0;
	double nmc;
	double fuvh_a;
	double fuvh_b;
	double fuvh_c;
	double irh_a;
	double irh_b;
	double crh;

	std::vector<double> cxhi_T;
	std::vector<double> cxhi_rate;
	std::vector<double> cxhi_rate2;
	double cxhi_minSlope, cxhi_maxSlope;

	std::vector<double> hr_T;
	std::vector<double> hr_rate;
	std::vector<double> hr_rate2;
	double hr_minSlope, hr_maxSlope;

	void initCollisionalExcitationHI(const Scalings& scale);
	void initRecombinationHII(const Scalings& scale);
	double collisionalExcitationHI(const double& nH, const double& HIIFRAC, const double& T);
	double recombinationHII(const double& nH, const double& HIIFRAC, const double& T);
	double coolingProcesses(const double& nH, const double& HIIFRAC, const double& T);
	double heatingProcesses(const double& nH, const double& HIIFRAC, const double& T, const double& Api, const double& F_FUV, const double& Av_FUV);
	double heating(const double& nH, const double& HIIFRAC, const double& T, const double& Api, const double& F_FUV, const double& tau);
};
