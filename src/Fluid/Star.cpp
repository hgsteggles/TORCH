#include "Star.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "MPI/MPI_Wrapper.hpp"
#include "Torch/Constants.hpp"
#include "GridCell.hpp"

void Star::initialise(std::shared_ptr<Constants> c, StarParameters sp, Location containing_core, const Vec3& delta_x) {
	consts = std::move(c);

	std::array<double, 3> mod;

	for (int i = 0; i < 3; ++i)
		mod[i] = (sp.faceSnap[i] && (sp.position[i] == 0)) ? 0.5 : 0;

	xc = std::array<double, 3>{ sp.position[0] + 0.5 - mod[0], sp.position[1] + 0.5 - mod[1], sp.position[2] + 0.5 - mod[2] };

	if (consts->nd < 3)
		xc[2] = 0;
	if (consts->nd < 2)
		xc[1] = 0;

	if (!sp.on)
		return;

	on = true;
	dx = delta_x;
	core = containing_core;
	photonEnergy = sp.photonEnergy;
	photonRate = sp.photonRate;
	windCellRadius = sp.windCellRadius;
	massLossRate = sp.massLossRate;
	windVelocity = sp.windVelocity;
	windTemperature = sp.windTemperature;
}

void Star::setWindCells(const CellContainer& cells) {
	double volume = 0;
	for (GridCell& cell : cells)
		volume += cell.vol;
	volume = MPIW::Instance().sum(volume);
	mdot = massLossRate/volume;
	edot = 0.5*mdot*windVelocity*windVelocity;
}

int Star::getWindCellRadius() const {
	return windCellRadius;
}

void Star::injectEnergyMomentum(CellContainer& windCells) {
	if (mdot != 0) {
		for (GridCell& cell : windCells) {
			//for (int idim = 0; idim < consts->nd; ++idim)
				//cell.UDOT[UID::VEL+idim] -= (mdot/cell.U[UID::DEN])*cell.U[UID::VEL+idim];
			cell.UDOT[UID::DEN] += mdot;
			cell.UDOT[UID::PRE] += edot;
			cell.UDOT[UID::HII] += (cell.U[UID::HII]/cell.U[UID::DEN])*mdot;
			cell.UDOT[UID::ADV] += mdot;
		}
	}
}

void Star::fixDensityPressure(CellContainer& windCells) {
	if (mdot != 0) {
		for (GridCell& cell : windCells) {
			//std::cout << "wind cell: xc[0] = " << cell.xc[0] << std::endl;
			double dist2 = 0;
			for (int idim = 0; idim < consts->nd; ++idim)
				dist2 += (cell.xc[idim] - xc[idim])*(cell.xc[idim] - xc[idim])*dx[idim]*dx[idim];
			if (dist2 == 0)
				dist2 = dx[0]*dx[0];
			double hii = cell.U[UID::HII]/cell.U[UID::DEN];

			cell.U[UID::DEN] = massLossRate/(4.0*consts->pi*dist2*windVelocity);
			//std::cout << cell.U[iden] << std::endl;
			cell.U[UID::HII] = hii*cell.U[UID::DEN];
			double ke = 0;
			for (int idim = 0; idim < consts->nd; ++idim) {
				cell.U[UID::VEL+idim] = windVelocity*cell.U[UID::DEN]*(cell.xc[idim] - xc[idim])*dx[idim]/std::sqrt(dist2);
				ke += 0.5*cell.U[UID::VEL+idim]*cell.U[UID::VEL+idim]/cell.U[UID::DEN];
			}
			double pre = (consts->specificGasConstant*(1.0 + cell.U[UID::HII]/cell.U[UID::DEN]))*cell.U[UID::DEN]*windTemperature;
			cell.U[UID::PRE] = pre/(cell.heatCapacityRatio-1) + ke;
			cell.U[UID::HII] = cell.U[UID::DEN];
		}
	}
}

