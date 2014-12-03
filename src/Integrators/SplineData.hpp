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

#include <vector>

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
	virtual ~SplineData() {};

	virtual double interpolate(double x) const = 0;
	static std::vector<double> spline(const std::vector<double>& x, const std::vector<double>& f, double yp1, double ypn);
	static double splint(const std::vector<double>& x, const std::vector<double>& f, const std::vector<double>& f2, double x2);
protected:
	std::vector<double> m_x;
	std::vector<double> m_y;
	std::vector<double> m_y2;
	double m_minSlope = 0;
	double m_maxSlope = 0;
};

class LinearSplineData : public SplineData {
public:
	LinearSplineData(const std::vector<std::pair<double, double>>& data);
	~LinearSplineData() {};

	virtual double interpolate(double x) const;
};

class LogSplineData : public SplineData {
public:
	LogSplineData(const std::vector<std::pair<double, double>>& data);
	~LogSplineData() {};

	virtual double interpolate(double x) const;
};



#endif /* SPLINEDATA_HPP_ */
