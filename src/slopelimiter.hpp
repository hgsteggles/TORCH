/**
 * 
 * @file slopelimiter.hpp
 *
 *  @author "Harrison Steggles"
 *  @date
 */

#ifndef SLOPELIMITER_HPP_
#define SLOPELIMITER_HPP_

namespace SlopeLimiter {
	double monotonizedCentral(double slopeA, double slopeB);
	double superbee(double slopeA, double slopeB);
	double minMod(double slopeA, double slopeB);
	double maxMod(double slopeA, double slopeB);
	double SF(double slopeA, double slopeB);
}

#endif /* SLOPELIMITER_HPP_ */
