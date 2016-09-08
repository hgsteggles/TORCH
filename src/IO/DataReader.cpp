/*
 * DataReader.cpp
 *
 *  Created on: 4 Mar 2015
 *      Author: harry
 */

#include <math.h>
#include <limits>
#include <algorithm>
#include <fstream>

#include "DataReader.hpp"
#include "Fluid/Fluid.hpp"
#include "Torch/Parameters.hpp"
#include "MPI/MPI_Wrapper.hpp"

DataParameters DataReader::readDataParameters(const std::string& filename) {
	DataParameters dp;
	MPIW& mpihandler = MPIW::Instance();

	mpihandler.serial([&] () {
		if (filename.compare("") != 0) {
			std::ifstream myfile(filename, std::ios_base::in);
			if (!myfile)
				throw std::runtime_error("DataReader::readDataParameters: invalid input file " + filename + ".");
			myfile >> dp.time >> dp.ncells[0] >> dp.ncells[1] >> dp.ncells[2];

			if (dp.ncells[1] == 1)
				dp.nd = 1;
			else if (dp.ncells[2] == 1)
				dp.nd = 2;
			else
				dp.nd = 3;

			std::array<double, 3> x0, x1;
			for (int i = 0; i < 3; ++i)
				x0[i] = x1[i] = 0;

			for (int idim = 0; idim < dp.nd; ++idim)
				myfile >> x0[idim];

			std::string ignore;

			myfile >> ignore >> ignore >> ignore;
			for (int idim = 0; idim < dp.nd; ++idim)
				myfile >> ignore;

			for (int idim = 0; idim < dp.nd; ++idim)
				myfile >> x1[idim];

			std::vector<double> dx;

			for (int i = 0; i < dp.nd; ++i) {
				if (std::abs(x1[i] - x0[i]) > std::numeric_limits<double>::min())
					dx.push_back(std::abs(x1[i] - x0[i]));
			}

			dp.dx = *(std::min_element(std::begin(dx), std::end(dx)));

			dp.sideLength = dp.dx*dp.ncells[0];

			myfile.close();
		}
	});

	return dp;
}

void DataReader::readGrid(const std::string& filename, const DataParameters& dp, Fluid& fluid) {
	MPIW& mpihandler = MPIW::Instance();
	Grid& grid = fluid.getGrid();

	mpihandler.serial([&] () {
		if (filename.compare("") != 0) {
			std::ifstream myfile(filename, std::ios_base::in);
			if (!myfile)
				throw std::runtime_error("DataReader::readGrid: invalid input file " + filename + ".");

			int nc = fluid.getGrid().ncells[0]*fluid.getGrid().ncells[1]*fluid.getGrid().ncells[2];

			std::string ignore;
			myfile >> ignore >> ignore >> ignore >> ignore;

			for (int i = 0; i < nc; ++i) {
				std::array<int, 3> xc;
				for (int idim = 0; idim < dp.nd; ++idim) {
					double x;
					myfile >> x;
					xc[idim] = (x/dp.dx);
				}
				for (int idim = dp.nd; idim < 3; ++idim)
					xc[idim] = 0;

				double density, pressure, hii;
				myfile >> density;
				myfile >> pressure;
				myfile >> hii;

				std::array<double, 3> velocity;
				for (int idim = 0; idim < dp.nd; ++idim)
					myfile >> velocity[idim];



				if (xc[0] < fluid.getGrid().getLeftX() || xc[0] > fluid.getGrid().getRightX())
					continue;

				int cellID = grid.locate((int)xc[0], (int)xc[1], (int)xc[2]);
				if (cellID != -1) {
					GridCell& cell = grid.getCell(cellID);

					cell.Q[UID::DEN] = density;
					cell.Q[UID::PRE] = pressure;
					cell.Q[UID::HII] = hii;

					for (int idim = 0; idim < dp.nd; ++idim)
						cell.Q[UID::VEL+idim] = velocity[idim];
					cell.heatCapacityRatio = fluid.heatCapacityRatio;
				}
			}
			myfile.close();
		}
	});
}

void DataReader::patchGrid(const std::string& filename, const std::array<int, 3>& offset, Fluid& fluid ) {
	MPIW& mpihandler = MPIW::Instance();
	Grid& grid = fluid.getGrid();

	DataParameters dp = readDataParameters(filename);

	mpihandler.serial([&] () {
		if (filename.compare("") != 0) {
			int nd = (fluid.getGrid().ncells[1] == 1) ? 1 : (fluid.getGrid().ncells[2] == 1) ? 2 : 3;
			if (dp.nd != nd)
				throw std::runtime_error("DataReader::patchGrid: patch ndims != grid ndims.");

			// How many patch cells fit into a grid cell?
			double ratio = fluid.getGrid().dx[0]/dp.dx;
			if (std::abs(ratio -(int)ratio) > std::numeric_limits<double>::min())
				throw std::runtime_error("DataReader::patchGrid: ratio of cell sizes not an integer.");

			for (int i = 0; i < 3; ++i) {
				double gridcellscovered = dp.ncells[i]/ratio;
				if (std::abs(gridcellscovered -(int)gridcellscovered) > std::numeric_limits<double>::min())
					throw std::runtime_error("DataReader::patchGrid: number of cells covered by patch not an integer.");
			}

			int r = std::pow((int)(ratio + 0.5), dp.nd);

			std::ifstream myfile(filename, std::ios_base::in);
			if (!myfile)
				throw std::runtime_error("DataReader::patchGrid: invalid input file " + filename + ".");

			std::string ignore;

			myfile >> ignore >> ignore >> ignore >> ignore;

			int npatchcells = dp.ncells[0]*dp.ncells[1]*dp.ncells[2];

			for (int i = 0; i < npatchcells; ++i) {
				std::array<int, 3> xc;
				for (int idim = 0; idim < dp.nd; ++idim) {
					double x;
					myfile >> x;
					xc[idim] = (x/fluid.getGrid().dx[0]) + offset[idim];
				}
				for (int idim = dp.nd; idim < 3; ++idim) {
					xc[idim] = 0;
				}
				myfile >> ignore >> ignore >> ignore;
				for (int idim = 0; idim < dp.nd; ++idim)
					myfile >> ignore;

				if (xc[0] < fluid.getGrid().getLeftX() || xc[0] > fluid.getGrid().getRightX())
					continue;

				int cellID = grid.locate((int)xc[0], (int)xc[1], (int)xc[2]);
				if (cellID != -1) {
					GridCell& cell = grid.getCell(cellID);
					cell.Q[UID::DEN] = 0;
					cell.Q[UID::PRE] = 0;
					cell.Q[UID::HII] = 0;
					for (int idim = 0; idim < dp.nd; ++idim)
						cell.Q[UID::VEL+idim] = 0;
				}
			}

			myfile.close();
			myfile.open(filename, std::ios_base::in);

			myfile >> ignore >> ignore >> ignore >> ignore;

			for (int i = 0; i < npatchcells; ++i) {
				std::array<int, 3> xc;
				for (int idim = 0; idim < dp.nd; ++idim) {
					double x;
					myfile >> x;
					xc[idim] = (x/fluid.getGrid().dx[0]) + offset[idim];
				}
				for (int idim = dp.nd; idim < 3; ++idim) {
					xc[idim] = 0;
				}

				double density, pressure, hii;
				myfile >> density;
				myfile >> pressure;
				myfile >> hii;

				std::array<double, 3> velocity;
				for (int idim = 0; idim < dp.nd; ++idim)
					myfile >> velocity[idim];

				if (xc[0] < fluid.getGrid().getLeftX() || xc[0] > fluid.getGrid().getRightX())
					continue;

				int cellID = grid.locate((int)xc[0], (int)xc[1], (int)xc[2]);
				if (cellID != -1) {
					GridCell& cell = grid.getCell(cellID);
					cell.Q[UID::DEN] += density/r;
					cell.Q[UID::PRE] += pressure/r;
					cell.Q[UID::HII] += hii/r;

					for (int idim = 0; idim < dp.nd; ++idim)
						cell.Q[UID::VEL+idim] += velocity[idim]/r;
				}
			}
			myfile.close();
		}
	});
}


