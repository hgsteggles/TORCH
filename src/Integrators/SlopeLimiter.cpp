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

double FalleLimiter::calculate(double a, double b) const {
	return(b*a < 1.0e-12 ? 0 : a*b*(a + b)/ (a*a + b*b));
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
	else if (type.compare("falle") == 0)
		return std::unique_ptr<SlopeLimiter>(new FalleLimiter());
	else if (type.compare("default") == 0)
		return std::unique_ptr<SlopeLimiter>(new FalleLimiter());
	else
		throw std::runtime_error("SlopeLimiterFactory::create: unknown type.");
}
