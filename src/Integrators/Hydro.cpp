/**
 * @file hydro.cpp
 */

#include "Hydro.hpp"
#include "Fluid.hpp"
//#include "grid3d.h"
#include "Grid.hpp"
#include "GridCell.hpp"
#include "Boundary.hpp"
#include "Star.hpp"
#include "MPI_Wrapper.hpp"
#include "Constants.hpp"

#include <string>
#include <limits> // numeric_limits
#include <cmath> // sqrt, etc
#include <iostream> // cout, cerr
#include <algorithm> //copy
#include <iterator> //begin, end

static void printQ(const FluidArray& Q) {
	std::cout << "density  = " << Q[UID::DEN] << std::endl;
	std::cout << "pressure = " << Q[UID::PRE] << std::endl;
	std::cout << "hii      = " << Q[UID::HII] << std::endl;
	std::cout << "vel0     = " << Q[UID::VEL+0] << std::endl;
	std::cout << "vel1     = " << Q[UID::VEL+1] << std::endl;
	std::cout << "vel2     = " << Q[UID::VEL+2] << std::endl;
}

Hydrodynamics::Hydrodynamics()
: Integrator("Hydrodynamics")
{ }

void Hydrodynamics::preTimeStepCalculations(Fluid& fluid) const {
	for (GridCell& cell : fluid.getGrid().getCells())
		cell.setSoundSpeed(fluid.calcSoundSpeed(cell.heatCapacityRatio, cell.Q[UID::PRE], cell.Q[UID::DEN]));
}

void Hydrodynamics::initialise(std::shared_ptr<Constants> c) {
	m_consts = std::move(c);
	m_riemannSolver = std::move(RiemannSolverFactory::create("default", m_consts->nd));
	m_slopeLimiter = std::move(SlopeLimiterFactory::create("default"));
}

void Hydrodynamics::integrate(double dt, Fluid& fluid) const {
	updateBoundaries(fluid);
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
	// Gridcells.
	for (GridCell& cell : fluid.getGrid().getCells()){
		for (int dim = 0; dim < m_consts->nd; ++dim) {
			GridCell* left = cell.ljoin[dim]->lcell;
			GridCell* right = cell.rjoin[dim]->rcell;
			piecewiseLinear(left->Q, cell.Q, right->Q, cell.QL[dim], cell.QR[dim]);
		}
	}
	// Left Boundaries.
	for (int dim = 0; dim < m_consts->nd; dim++) {
		for (GridCell* cptr : fluid.getGrid().getLeftBoundaries()[dim]->getGhostCells()) {
			GridCell* left = cptr->left[dim];
			GridCell* right = cptr->rjoin[dim]->rcell;
			piecewiseLinear(left->Q, cptr->Q, right->Q, cptr->QL[dim], cptr->QR[dim]);
		}
	}
	// Right Boundaries.
	for (int dim = 0; dim < m_consts->nd; dim++) {
		for (GridCell* cptr : fluid.getGrid().getRightBoundaries()[dim]->getGhostCells()) {
			GridCell* left = cptr->ljoin[dim]->lcell;
			GridCell* right = cptr->right[dim];
			piecewiseLinear(left->Q, cptr->Q, right->Q, cptr->QL[dim], cptr->QR[dim]);
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
	for (GridCell& cell : fluid.getGrid().getCells()){
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
	if (fluid.getGrid().spatialOrder == 1)
		reconstruct(fluid);
	for (int dim = 0; dim < m_consts->nd; ++dim){
		for (GridJoin& join  : fluid.getGrid().getJoins(dim)) {
			if (fluid.getGrid().spatialOrder == 1) {
				//characteristicWaveSpeeds(S_l, S_r, join.lcell->QR[dim], join.rcell->QL[dim], join.lcell->heatCapacityRatio, dim);
				double a_l2 = soundSpeedSqrd(join.lcell->QR[dim][UID::PRE], join.lcell->QR[dim][UID::DEN], join.lcell->heatCapacityRatio);
				double a_r2 = soundSpeedSqrd(join.rcell->QL[dim][UID::PRE], join.rcell->QL[dim][UID::DEN], join.rcell->heatCapacityRatio);
				m_riemannSolver->solve(join.F, join.lcell->QR[dim], join.rcell->QL[dim], a_l2, a_r2, join.lcell->heatCapacityRatio, dim);
			}
			else if (fluid.getGrid().spatialOrder == 0) {
				//characteristicWaveSpeeds(S_l, S_r, join.lcell->Q, join.rcell->Q, join.lcell->heatCapacityRatio, dim);
				double a_l2 = soundSpeedSqrd(join.lcell->Q[UID::PRE], join.lcell->Q[UID::DEN], join.lcell->heatCapacityRatio);
				double a_r2 = soundSpeedSqrd(join.rcell->Q[UID::PRE], join.rcell->Q[UID::DEN], join.rcell->heatCapacityRatio);
				m_riemannSolver->solve(join.F, join.lcell->Q, join.rcell->Q, a_l2, a_r2, join.lcell->heatCapacityRatio, dim);
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

	for (GridCell& cell : fluid.getGrid().getCells()) {
		//Fluxes.
		for (int dim = 0; dim < m_consts->nd; ++dim) {
			for (int i = 0; i < UID::N; ++i) {
				if (cell.ljoin[dim] == nullptr || cell.rjoin[dim] == nullptr)
					throw std::runtime_error("Hydrodynamics::updateSourceTerms: ljoin[dim] or rjoin[dim] is nullptr." + cell.printInfo());
				cell.UDOT[i] += (cell.ljoin[dim]->area/cell.vol)*cell.ljoin[dim]->F[i];
				cell.UDOT[i] -= (cell.rjoin[dim]->area/cell.vol)*cell.rjoin[dim]->F[i];
			}
		}

		//Internal Energy Source.
		double div_v = 0;
		if (fluid.getGrid().spatialOrder == 1) {
			for (int i = 0; i < m_consts->nd; ++i)
				div_v += (cell.QR[i][UID::VEL+i]-cell.QL[i][UID::VEL+i])/fluid.getGrid().dx[i];
		}
		else {
			for (int i = 0; i < m_consts->nd; ++i)
				div_v += m_slopeLimiter->calculate(cell.Q[UID::VEL+i] - cell.ljoin[i]->lcell->Q[UID::VEL+i], cell.rjoin[i]->rcell->Q[UID::VEL+i] - cell.Q[UID::VEL+i])/fluid.getGrid().dx[i];
		}

		//Geometric.
		if (rzp) {
			double r = fluid.getGrid().dx[0]*cell.xc[0];
			cell.UDOT[UID::VEL+0] += cell.Q[UID::PRE]/r;
			div_v += cell.Q[UID::VEL+0]/r;
		}
		else if (rtp) {
			double r, area1, area2;
			area1 = cell.rjoin[0]->area;
			area2 = cell.ljoin[0]->area;
			r = cell.vol/(area1 - area2);
			cell.UDOT[UID::VEL+0] += cell.Q[UID::PRE]/r;
			div_v += 2.0*cell.Q[UID::VEL+0]/r;
		}

		//Internal Energy Source.
		//cell.UDOT[UID::EINT] -= cell.Q[UID::PRE]*div_v;
	}

	if (fluid.getStar().on)
		fluid.getStar().injectEnergyMomentum((fluid.getGrid().getWindCells()));
}

void Hydrodynamics::updateBoundaries(Fluid& fluid) const {
	if (MPIW::Instance().getRank()%2 == 0) {
		fluid.getGrid().getLeftBoundaries()[0]->applyBC();
		fluid.getGrid().getRightBoundaries()[0]->applyBC();
	}
	else {
		fluid.getGrid().getRightBoundaries()[0]->applyBC();
		fluid.getGrid().getLeftBoundaries()[0]->applyBC();
	}
	for (int i = 1; i < m_consts->nd; ++i) {
		fluid.getGrid().getRightBoundaries()[i]->applyBC();
		fluid.getGrid().getLeftBoundaries()[i]->applyBC();
	}
	/*
	mpihandler.serial([&] () {
		std::cout << std::endl;
		if (MPIW::Instance().getRank() == 0) {
			GridCell& cell0 = **fluid.getGrid().getRightBoundaries()[0]->getGhostCells().begin();
			GridCell& cell1 = *cell0.right[0];
			GridCell& c0 = *cell0.left[0];
			GridCell& c1 = *cell1.left[0];
			std::cout << cell0.printCoords() << cell0.Q[UID::PRE] << std::endl;
			std::cout << cell1.printCoords() << cell1.Q[UID::PRE] << std::endl;
			std::cout << c0.printCoords() << c0.Q[UID::PRE] << std::endl;
			std::cout << c1.printCoords() << c1.Q[UID::PRE] << std::endl;
		}
		else if (MPIW::Instance().getRank() == 1) {
			GridCell& cell0 = **fluid.getGrid().getLeftBoundaries()[0]->getGhostCells().begin();
			GridCell& cell1 = *cell0.left[0];
			GridCell& c0 = *cell0.right[0];
			GridCell& c1 = *cell1.right[0];
			std::cout << cell0.printCoords() << cell0.Q[UID::PRE] << std::endl;
			std::cout << cell1.printCoords() << cell1.Q[UID::PRE] << std::endl;
			std::cout << c0.printCoords() << c0.Q[UID::PRE] << std::endl;
			std::cout << c1.printCoords() << c1.Q[UID::PRE] << std::endl;
		}
		std::cout << std::endl;
	});
	*/
}
