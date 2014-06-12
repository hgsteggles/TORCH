/**
 * Provides the Thermodynamics class.
 * @file thermodynamics.hpp
 * @author Harrison Steggles
 * @date 12/06/2014 - The first version.
 */

#include <vector>
#include <cmath>

/**
 * @class Thermodynamics
 * @brief Contains parameters and methods for calculating heating and cooling due to atomic processes in the gas around an ionizing star.
 * @version 0.1, 12/06/2014
 */
class Thermodynamics {
public:
	Thermodynamics();
	double fluxFUV(const double& Q_FUV, const double dist_sqrd);
	double heating(const double& nH, const double& HIIFRAC, const double& T,  const double& F_FUV, const double& tau);
private:
	double sigmaV; //<! Dust extinction cross-section [cm2 H-1].

	std::vector<double> cxhi_T;
	std::vector<double> cxhi_rate;
	std::vector<double> cxhi_rate2;
	double cxhi_minSlope, cxhi_maxSlope;

	std::vector<double> hr_T;
	std::vector<double> hr_rate;
	std::vector<double> hr_rate2;
	double hr_minSlope, hr_maxSlope;

	void initCollisionalExcitationHI();
	void initRecombinationHII();
	double collisionalExcitationHI(const double& T);
	double recombinationHII(const double& T);
	double coolingProcesses(const double& nH, const double& HIIFRAC, const double& T);
	double heatingProcesses(const double& nH, const double& T,  const double& F_FUV, const double& Av_FUV);
};

namespace Cooling {



}
