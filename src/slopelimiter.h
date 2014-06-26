/**
 * Provides methods for calculating the gradient of the fluid variables in a GridCell.
 *
 * @file star.h
 *
 * @author Harrison Steggles
 *
 * @date 12/06/2014 - The first version.
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
