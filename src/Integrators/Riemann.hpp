/** Provides the RiemannSolver base abstract class and HartenLaxLeerSolver, HartenLaxLeerContactSolver and
 * RotatedHartenLaxLeerSolver subclasses.
 *
 * @file Riemann.hpp
 *
 * @author Harrison Steggles
 */

#ifndef RIEMANN_HPP_
#define RIEMANN_HPP_

#include <memory>
#include <string>

#include "Torch/Common.hpp"

/**
 * @class RiemannSolver
 *
 * @brief A base class for solving the Riemann problem.
 */
class RiemannSolver {
public:
	RiemannSolver(int nd);
	virtual ~RiemannSolver() { };
	virtual void solve(FluidArray& F, const FluidArray& Q_l, const FluidArray& Q_r, double a_l2, double a_r2, double gamma, int dim) const = 0;
	int getNumberDimensions() const;
private:
	int m_ND;
};

/**
 * @class HartenLaxLeerContactSolver
 *
 * @brief RiemannSolver subclass. Refer to "Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction".
 */
class HartenLaxLeerContactSolver : public RiemannSolver {
public:
	HartenLaxLeerContactSolver(int nd);
	virtual void solve(FluidArray& F, const FluidArray& Q_l, const FluidArray& Q_r, double a_l2, double a_r2, double gamma, int dim) const;
};

/**
 * @class HartenLaxLeerSolver
 *
 * @brief RiemannSolver subclass. Refer to "Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction".
 */
class HartenLaxLeerSolver : public RiemannSolver {
public:
	HartenLaxLeerSolver(int nd);
	virtual void solve(FluidArray& F, const FluidArray& Q_l, const FluidArray& Q_r, double a_l2, double a_r2, double gamma, int dim) const;
};

/**
 * @class RotatedHartenLaxLeerSolver
 *
 * @brief RiemannSolver subclass. Refer to "Nishikawa, H. & Kitamura, K. 2008, Journal of Computational Physics, 227, 2560".
 */
class RotatedHartenLaxLeerSolver : public RiemannSolver {
public:
	RotatedHartenLaxLeerSolver(int nd);
	virtual void solve(FluidArray& F, const FluidArray& Q_l, const FluidArray& Q_r, double a_l2, double a_r2, double gamma, int dim) const;
private:
	HartenLaxLeerContactSolver m_hllc;
	HartenLaxLeerSolver m_hll;
};

/**
 * @class RiemannSolverFactory
 *
 * @brief Creates RiemannSolvers.
 */
class RiemannSolverFactory {
public:
	static std::unique_ptr<RiemannSolver> create(const std::string& type, int ndims);
};


#endif // RIEMANN_HPP_
