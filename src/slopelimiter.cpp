/**
 * @file slopelimiter.cpp
 */

#include "slopelimiter.h"

#include <cmath>

double SlopeLimiter::monotonizedCentral(double slopeA, double slopeB) {
	return minMod((slopeA+slopeB)/2.0, minMod(2.0*slopeA, 2.0*slopeB));
}

double SlopeLimiter::superbee(double slopeA, double slopeB) {
	return maxMod(minMod(slopeB, 2.0*slopeA), minMod(2.0*slopeB, slopeA));
}

double SlopeLimiter::minMod(double slopeA, double slopeB) {
	if (slopeA*slopeB > 0)
		return std::abs(slopeA) < std::abs(slopeB) ? slopeA : slopeB;
	else
		return 0;
}

double SlopeLimiter::maxMod(double slopeA, double slopeB) {
	if (slopeA*slopeB > 0)
		return std::abs(slopeA) > std::abs(slopeB) ? slopeA : slopeB;
	else
		return 0;
}


double SlopeLimiter::SF(double slopeA, double slopeB) {
  return(slopeB*slopeA < 1.0e-12 ? 0 : slopeA*slopeB*(slopeA + slopeB)/ (slopeA*slopeA + slopeB*slopeB));
}


