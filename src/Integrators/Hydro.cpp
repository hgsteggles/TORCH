#include "Hydro.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

#include "Fluid/Fluid.hpp"
#include "Fluid/Grid.hpp"
#include "Fluid/GridCell.hpp"
#include "Fluid/Star.hpp"
#include "MPI/MPI_Wrapper.hpp"
#include "Torch/Constants.hpp"

Hydrodynamics::Hydrodynamics()
: Integrator("Hydrodynamics")
{ }

void Hydrodynamics::preTimeStepCalculations(Fluid& fluid) const {
	for (GridCell& cell : fluid.getGrid().getIterable("GridCells"))
		cell.setSoundSpeed(fluid.calcSoundSpeed(cell.heatCapacityRatio, cell.Q[UID::PRE], cell.Q[UID::DEN]));
}

void Hydrodynamics::initialise(std::shared_ptr<Constants> c) {
	m_consts = std::move(c);
	m_riemannSolver = std::move(RiemannSolverFactory::create("default", m_consts->nd));
	m_slopeLimiter = std::move(SlopeLimiterFactory::create("default"));
}

void Hydrodynamics::integrate(double dt, Fluid& fluid) const {
	fluid.getGrid().applyBCs();
	calcFluxes(fluid);
}

void Hydrodynamics::setRiemannSolver(std::unique_ptr<RiemannSolver> riemannSolver) {
	m_riemannSolver = std::move(riemannSolver);
}

void Hydrodynamics::setSlopeLimiter(std::unique_ptr<SlopeLimiter> slopeLimiter) {
	m_slopeLimiter = std::move(slopeLimiter);
}

void Hydrodynamics::piecewiseLinear(FluidArray& Q_l, FluidArray& Q_c, FluidArray& Q_r, FluidArray& left_interp, FluidArray& right_interp) const {
	for (int iq = 0; iq < UID::N; ++iq) {
		double dQdr = m_slopeLimiter->calculate(Q_c[iq] - Q_l[iq], Q_r[iq] - Q_c[iq]);
		left_interp[iq] = Q_c[iq] - 0.5*dQdr;
		right_interp[iq] = Q_c[iq] + 0.5*dQdr;
	}
}

void Hydrodynamics::reconstruct(Fluid& fluid) const {
	for (GridCell& cell : fluid.getGrid().getIterable("GridCells")) {
		for (int dim = 0; dim < m_consts->nd; ++dim) {
			GridCell& left = fluid.getGrid().left(dim, cell);
			GridCell& right = fluid.getGrid().right(dim, cell);
			piecewiseLinear(left.Q, cell.Q, right.Q, cell.QL[dim], cell.QR[dim]);
		}
	}
	for (int dim = 0; dim < m_consts->nd; ++dim) {
		for (GridCell& cell : fluid.getGrid().getIterable("GhostCells"+std::to_string(dim))) {
			GridCell& left = fluid.getGrid().left(dim, cell);
			GridCell& right = fluid.getGrid().right(dim, cell);
			piecewiseLinear(left.Q, cell.Q, right.Q, cell.QL[dim], cell.QR[dim]);
		}
	}
}

double Hydrodynamics::soundSpeedSqrd(const double pre, const double den, const double gamma) const {
	return gamma*pre/den;
}

double Hydrodynamics::soundSpeed(const double pre, const double den, const double gamma) const {
	return std::sqrt(soundSpeedSqrd(pre, den, gamma));
}

double Hydrodynamics::calculateTimeStep(double dt_max, Fluid& fluid) const {
	double tmin = dt_max;
	for (GridCell& cell : fluid.getGrid().getIterable("GridCells")) {
		double inv_t = 0;
		double ss = soundSpeed(cell.Q[UID::PRE], cell.Q[UID::DEN], fluid.heatCapacityRatio);
		for (int dim = 0; dim < m_consts->nd; ++dim)
			inv_t += (fabs(cell.Q[UID::VEL+dim]) + ss)/fluid.getGrid().dx[dim];
		if (inv_t == 0)
			inv_t = 1.0/dt_max;
		tmin = std::min(tmin, 0.5/inv_t);
	}
	return std::min(tmin, dt_max);
}

void Hydrodynamics::calcFluxes(Fluid& fluid) const {
	Grid& grid = fluid.getGrid();

	if (grid.spatialOrder == 1)
		reconstruct(fluid);
	for (int dim = 0; dim < m_consts->nd; ++dim){
		for (GridJoin& join  : grid.getJoins(dim)) {
			GridCell& left = grid.joinedLeft(join);
			GridCell& right = grid.joinedRight(join);
			if (fluid.getGrid().spatialOrder == 1) {
				if (left.QR[dim][UID::DEN] <= 1e-20 || right.QL[dim][UID::DEN] <= 1e-20) {
					std::cout << "left: " << left.xc[0] << " " << left.xc[1] << " " << left.xc[2] << std::endl;
					std::cout << "right: " << right.xc[0] << " " << right.xc[1] << " " << right.xc[2] << std::endl;
				}
				double a_l2 = soundSpeedSqrd(left.QR[dim][UID::PRE], left.QR[dim][UID::DEN], left.heatCapacityRatio);
				double a_r2 = soundSpeedSqrd(right.QL[dim][UID::PRE], right.QL[dim][UID::DEN], right.heatCapacityRatio);
				m_riemannSolver->solve(join.F, left.QR[dim], right.QL[dim], a_l2, a_r2, left.heatCapacityRatio, dim);
			}
			else if (fluid.getGrid().spatialOrder == 0) {
				double a_l2 = soundSpeedSqrd(left.Q[UID::PRE], left.Q[UID::DEN], left.heatCapacityRatio);
				double a_r2 = soundSpeedSqrd(right.Q[UID::PRE], right.Q[UID::DEN], right.heatCapacityRatio);
				m_riemannSolver->solve(join.F, left.Q, right.Q, a_l2, a_r2, left.heatCapacityRatio, dim);
			}
			else {
				throw std::runtime_error("Hydrodynamics::calcFluxes: invalid order(=" + std::to_string(fluid.getGrid().spatialOrder) + "). Valid orders = {0, 1}.");
			}
		}
	}
}

void Hydrodynamics::updateSourceTerms(double dt, Fluid& fluid) const {
	bool rtp = (fluid.getGrid().geometry == Geometry::SPHERICAL);
	bool rzp = (fluid.getGrid().geometry == Geometry::CYLINDRICAL);

	Grid& grid = fluid.getGrid();

	for (GridCell& cell : grid.getIterable("GridCells")) {
		//Fluxes.
		for (int dim = 0; dim < m_consts->nd; ++dim) {
			for (int i = 0; i < UID::N; ++i) {
				if (!grid.joinExists(dim, cell.ljoinID) || !grid.joinExists(dim, cell.rjoinID))
					throw std::runtime_error("Hydrodynamics::updateSourceTerms: ljoin[dim] or rjoin[dim] is nullptr." + cell.printInfo());
				cell.UDOT[i] += (grid.leftJoin(dim, cell).area/cell.vol)*grid.leftJoin(dim, cell).F[i];
				cell.UDOT[i] -= (grid.rightJoin(dim, cell).area/cell.vol)*grid.rightJoin(dim, cell).F[i];
			}
		}

		for (int i = 0; i < 3; ++i) {
			cell.UDOT[UID::VEL+i] += cell.GRAV[i];
			cell.UDOT[UID::PRE] += cell.Q[UID::VEL+i]*cell.GRAV[i];
		}

		//Geometric.
		if (rzp) {
			double r = fluid.getGrid().dx[0]*cell.xc[0];
			cell.UDOT[UID::VEL+0] += cell.Q[UID::PRE]/r;
		}
		else if (rtp) {
			double r, area1, area2;
			area1 = grid.rightJoin(0, cell).area;
			area2 = grid.leftJoin(0, cell).area;
			r = cell.vol/(area1 - area2);
			cell.UDOT[UID::VEL+0] += cell.Q[UID::PRE]/r;
		}

	}

	if (fluid.getStar().on)
		fluid.getStar().injectEnergyMomentum(grid);
}
