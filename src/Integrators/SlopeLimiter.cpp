#include "SlopeLimiter.hpp"

#include <cmath>
#include <stdexcept>

inline double minmod(double a, double b) {
	if (a*b > 0)
		return std::abs(a) < std::abs(b) ? a : b;
	else
		return 0;
}

inline double maxmod(double a, double b) {
	if (a*b > 0)
		return std::abs(a) > std::abs(b) ? a : b;
	else
		return 0;
}

double MinmodLimiter::calculate(double a, double b) const {
	return minmod(a, b);
}

double MaxmodLimiter::calculate(double a, double b) const {
	return maxmod(a, b);
}

double MonotonisedCentralLimiter::calculate(double a, double b) const {
	return minmod((a+b)/2.0, minmod(2.0*a, 2.0*b));
}

double SuperbeeLimiter::calculate(double a, double b) const {
	return maxmod(minmod(b, 2.0*a), minmod(2.0*b, a));
}

double LeerLimiter::calculate(double a, double b) const {
	if (b == 0)
		return 0;
	else {
		double r = std::abs(a/b);
		return (a + b*r)/(1.0 + r);
	}
}

double OspreLimiter::calculate(double a, double b) const {
	if (a*b <= 1.0e-30)
		return 0;
	else
		return 1.5*b*(a*a + a*b)/(a*a + a*b + b*b);
}

double AlbadaLimiter::calculate(double a, double b) const {
	return a*b <= 1.0e-30 ? 0 : a*b*(a + b) / (a*a + b*b);
}

std::unique_ptr<SlopeLimiter> SlopeLimiterFactory::create(std::string type) {
	if (type.compare("monotonised_central") == 0)
		return std::unique_ptr<SlopeLimiter>(new MonotonisedCentralLimiter());
	else if (type.compare("superbee") == 0)
		return std::unique_ptr<SlopeLimiter>(new SuperbeeLimiter());
	else if (type.compare("minmod") == 0)
		return std::unique_ptr<SlopeLimiter>(new MinmodLimiter());
	else if (type.compare("maxmod") == 0)
		return std::unique_ptr<SlopeLimiter>(new MaxmodLimiter());
	else if (type.compare("leer") == 0)
		return std::unique_ptr<SlopeLimiter>(new LeerLimiter());
	else if (type.compare("ospre") == 0)
		return std::unique_ptr<SlopeLimiter>(new OspreLimiter());
	else if (type.compare("albada") == 0)
		return std::unique_ptr<SlopeLimiter>(new AlbadaLimiter());
	else if (type.compare("default") == 0)
		return std::unique_ptr<SlopeLimiter>(new AlbadaLimiter());
	else
		throw std::runtime_error("SlopeLimiterFactory::create: unknown type.");
}
