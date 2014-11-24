/** Provides the SplineData class.
 *
 * @file SplineData.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef SPLINEDATA_HPP_
#define SPLINEDATA_HPP_

/**
 * @class SplineData
 *
 * @brief Contains cubic spline data for interpolation.
 *
 * Cubic spline interpolation and powerlaw extrapolation (using logarithm gradients at the extremes).
 *
 * @version 0.8, 24/11/2014
 */
class SplineData {
public:
	SplineData(const std::vector<std::pair<double, double>>& data);

	double interpolate(double x) const;
private:
	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<double> m_y2;
	double m_minSlope = 0;
	double m_maxSlope = 0;

	std::vector<double> spline(const std::vector<double>& x, const std::vector<double>& f, double yp1, double ypn);
	double splint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, double x2);
};



#endif /* SPLINEDATA_HPP_ */
