/** Provides the SlopeLimiter abstract base class and a few subclasses.
 *
 * @file SlopeLimiter.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef SLOPELIMITER_HPP_
#define SLOPELIMITER_HPP_

#include <functional>
#include <memory>
#include <iostream>

enum class LimiterName : unsigned int {MONOTONISED_CENTRAL, SUPERBEE, MINMOD, MAXMOD, FALLE};

/**
 * @class SlopeLimiter
 *
 * @brief A base class for calculating a slope which keeps the calling scheme total variation diminishing.
 *
 * @version 0.8, 24/11/2014
 */
class SlopeLimiter {
public:
	virtual ~SlopeLimiter() { };
	virtual double calculate(double a, double b) const = 0;
};

class MonotonisedCentralLimiter : public SlopeLimiter {
public:
	double calculate(double a, double b) const;
};

class SuperbeeLimiter : public SlopeLimiter {
public:
	double calculate(double a, double b) const;
};

class MinmodLimiter : public SlopeLimiter {
public:
	double calculate(double a, double b) const;
};

class MaxmodLimiter : public SlopeLimiter {
public:
	double calculate(double a, double b) const;
};

class LeerLimiter : public SlopeLimiter {
public:
	double calculate(double a, double b) const;
};

class OspreLimiter : public SlopeLimiter {
public:
	double calculate(double a, double b) const;
};

class AlbadaLimiter : public SlopeLimiter {
public:
	double calculate(double a, double b) const;
};

class SlopeLimiterFactory {
public:
	static std::unique_ptr<SlopeLimiter> create(std::string type);
};

#endif /* SLOPELIMITER_HPP_ */
