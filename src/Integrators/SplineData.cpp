#include "SplineData.hpp"

#include <cmath>
#include <stdexcept>

SplineData::SplineData(const std::vector<std::pair<double, double>>& data) {
	int ndata = 0;
	for (const std::pair<double, double>& dataPair : data) {
		for (double& xi : m_x) {
			if (xi == dataPair.first)
				throw std::runtime_error("SplineData::Constructor: x values must be distinct."); //The xa's must be distinct.
		}
		m_x.push_back(dataPair.first);
		m_y.push_back(dataPair.second);
		++ndata;
	}
	if (ndata == 0)
		throw std::runtime_error("SplineData::Constructor: zero length data vector passed to constructor.");

	m_y2 = spline(m_x, m_y, 1.e99, 1.e99);
	// Logarithmic slopes at either end of the domain.
	m_minSlope = (std::log10(m_y[1])-std::log10(m_y[0]))/(std::log10(m_x[1])-std::log10(m_x[0]));
	m_maxSlope = (std::log10(m_y[ndata-1])-std::log10(m_y[ndata-2]))/(std::log10(m_x[ndata-1])-std::log10(m_x[ndata-2]));
}

/**
 * @brief Calculates a cubic spline interpolated value from the data this object was initialised with.
 * @param x Interpolation location.
 * @return Cubic spline interpolated value.
 */
double SplineData::interpolate(double x) const {
	double rate = 0;
	if (x > m_x[m_x.size()-1])
		rate = m_y[m_x.size()-1] * std::pow( x/m_x[m_x.size()-1], m_maxSlope );
	else if ( x < m_x[0])
		rate = m_y[0] * std::pow(x/m_x[0], m_minSlope);
	else {
		int klo, khi, k;
		double h,b,a;
		int n = m_x.size();
		// We will find the right place in the table by means of bisection. This is
		// optimal if sequential calls to this routine are at random values of x.
		// If sequential calls are in order, and closely spaced, one would do better
		// to store previous values of klo and khi and test if they remain appropriate
		// on the next call.
		klo = 0;
		khi = n-1;

		while (khi-klo > 1) {
			k = (khi+klo) >> 1;
			if (m_x[k] > x)
				khi = k;
			else
				klo = k;
		}

		// klo and khi now bracket the input value of x.
		h = m_x[khi]-m_x[klo];

		a = (m_x[khi]-x)/h;
		b = (x-m_x[klo])/h; // Cubic spline polynomial is now evaluated.

		rate = a*m_y[klo]+b*m_y[khi]+((a*a*a-a)*m_y2[klo]+(b*b*b-b)*m_y2[khi])*(h*h)/6.0;
	}
	return rate;
}

/**
 * @brief Calculates second derivatives of an interpolating function.
 * @param x Vector of function arguments.
 * @param f Vector of function values given a vector of function arguments.
 * @param yp1 First derivative of function at x[0].
 * @param ypn First derivative of function at x[n-1].
 * @return Vector of second derivatives of interpolating function.
 */
std::vector<double> SplineData::spline(const std::vector<double>& x, const std::vector<double>& f, double yp1, double ypn){
	int n = x.size();
	std::vector<double> y2(n, 0.0);
	std::vector<double> u(n-1, 0.0);
	int i,k;
	double p,qn,sig,un;
	//The lower boundary condition is set either to be "natural"
	//or else to have a specified first derivative.
	if ( yp1 > 0.99e30) {
		y2[0] = 0.0;
		u[0] = 0.0;
	}
	else {
		y2[0] = -0.5;
		u[0] = (3.0/(x[1]-x[0]))*((f[1]-f[0])/(x[1]-x[0])-yp1);
	}
	//This is the decomposition loop of the tridiagonal algorithm.
	//y2 and u are used for temporary storage of the decomposed factors.
	for ( i = 1; i < n-1; i++ ) {
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sig*y2[i-1]+2.0;
		y2[i] = (sig-1.0)/p;
		u[i] = (f[i+1]-f[i])/(x[i+1]-x[i]) - (f[i]-f[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	//The upper boundary condition is set either to be "natural"
	// or else to have a specified first derivative.
	if ( ypn > 0.99e30 ) {
		qn = 0.0;
		un = 0.0;
	}
	else {
		qn = 0.5;
		un = (3.0/(x[n-1]-x[n-2]))*(ypn-(f[n-1]-f[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	//This is the backsubstitution loop of the tridiagonal algorithm.
	for ( k = n-2; k >= 0; k-- )
		y2[k] = y2[k]*y2[k+1]+u[k];
	return y2;
}
