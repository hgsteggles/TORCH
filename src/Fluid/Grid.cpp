#include "Grid.hpp"

#include "IO/ProgressBar.hpp"
#include "IO/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>


static int safe_round(double val) {
	if (val < 0)
		return (int)(val-0.5);
	else
		return (int)(val+0.5);
}

Grid::Grid()
{ }

GridCell& Grid::getCell(int id) {
	if (id < 0 || id >= (int)m_cells.size()) {
		throw std::runtime_error("Grid::getCell: out of bounds.\n");
	}
	return m_cells[id];
}

const GridCell& Grid::getCell(int id) const {
	if (id < 0 || id >= (int)m_cells.size()) {
		throw std::runtime_error("Grid::getCell: out of bounds.\n");
	}
	return m_cells[id];
}

Looper Grid::getIterable(const std::string& name) {
	return m_cellCollection.getIterable(name);
}

ConstLooper Grid::getIterable(const std::string& name) const {
	return m_cellCollection.getIterable(name);
}

int Grid::getLeftX() const {
	return m_leftX;
}

int Grid::getRightX() const {
	return m_rightX;
}

void Grid::setLeftX(int x) {
	m_leftX = x;
}

void Grid::setRightX(int x) {
	m_rightX = x;
}

std::vector<GridCell>& Grid::getCells() {
	return m_cells;
}

std::vector<GridJoin>& Grid::getJoins(int dim) {
	return m_joins[dim];
}

const std::vector<GridCell>& Grid::getCells() const {
	return m_cells;
}

const std::vector<GridJoin>& Grid::getJoins(int dim) const {
	return m_joins[dim];
}

PartitionManager& Grid::getPartitionManager() {
	return partition;
}

std::vector<int>& Grid::getCausalIndices() {
	return m_causalIndices;
}

std::vector<int>& Grid::getOrderedIndices(const std::string& name) {
	if (orderedIndices.find(name) == orderedIndices.end())
		Logger::Instance().print<SeverityType::WARNING>(name + " not found...\n");
	return (std::vector<int>&)orderedIndices[name];
}

std::vector<Bound>& Grid::getBoundaries() {
	return m_boundaries;
}

void Grid::addOrderedIndex(const std::string& name, int index) {
	std::vector<int>& order = (std::vector<int>&)orderedIndices[name];
	order.push_back(index);
}

bool Grid::cellExists(int id) const {
	return id >= 0 && id < (int)m_cells.size();
}

int Grid::left(int dim, int fromCellID) {
	return m_cells[fromCellID].leftID[dim];
}

int Grid::right(int dim, int fromCellID) {
	return m_cells[fromCellID].rightID[dim];
}

int Grid::ghostLeft(int dim, int fromCellID) {
	int joinID = m_cells[fromCellID].ljoinID[dim];
	if (joinID != -1)
		return m_joins[dim][joinID].lcellID;
	else
		return -1;
}

int Grid::ghostRight(int dim, int fromCellID) {
	return m_joins[dim][m_cells[fromCellID].rjoinID[dim]].rcellID;
}

GridCell& Grid::left(int dim, GridCell& fromCell) {
	return m_cells[fromCell.leftID[dim]];
}

GridCell& Grid::right(int dim, GridCell& fromCell) {
	return m_cells[fromCell.rightID[dim]];
}

GridCell& Grid::ghostLeft(int dim, GridCell& fromCell) {
	return m_cells[m_joins[dim][fromCell.ljoinID[dim]].lcellID];
}

GridCell& Grid::ghostRight(int dim, GridCell& fromCell) {
	return m_cells[m_joins[dim][fromCell.rjoinID[dim]].rcellID];
}

GridCell& Grid::joinedLeft(GridJoin& join) {
	return m_cells[join.lcellID];
}

GridCell& Grid::joinedRight(GridJoin& join) {
	return m_cells[join.rcellID];
}

GridJoin& Grid::leftJoin(int dim, GridCell& fromCell) {
	return m_joins[dim][fromCell.ljoinID[dim]];
}

GridJoin& Grid::rightJoin(int dim, GridCell& fromCell) {
	return m_joins[dim][fromCell.rjoinID[dim]];
}

int Grid::flatIndex(int ci, int cj, int ck) {
	return ci + coreCells[0]*(cj + coreCells[1]*ck);
}

Coords Grid::unflatCoords(int flat_index) {
	Coords coords;
	coords[2] = (int)(0.5 + flat_index/(coreCells[0]*coreCells[1]));
	coords[1] = (int)(0.5 + (flat_index - coords[2]*coreCells[0]*coreCells[1])/coreCells[0]);
	coords[0] = flat_index - coreCells[0]*(coords[1] + coreCells[1]*coords[2]);
	return coords;
}

int Grid::traverse1D(int dim, const int dxc, int fromCellID) {
	int id = fromCellID;
	int adxc = std::abs(dxc);
	for (int i = 0; i < adxc && id != -1; ++i) {
		id = dxc > 0 ? m_cells[id].rightID[dim] : m_cells[id].leftID[dim];
	}
	return id;
}

int Grid::traverse3D(const int d1, const int d2, const int d3, const int dc1, const int dc2, const int dc3, int fromCellID) {
	return traverse1D(d3, dc3, traverse1D(d2, dc2, traverse1D(d1, dc1, fromCellID)));
}

int Grid::traverseOverJoins1D(const int dim, const int dxc, int fromCellID) {
	int id = fromCellID;
	for(int i = 0; i < abs(dxc) && id != -1; ++i){
		if (dxc > 0 && m_cells[id].rjoinID[dim] != -1)
			id = ghostRight(dim, id);
		else if (dxc < 0 && m_cells[id].ljoinID[dim] != 0)
			id = ghostLeft(dim, id);
		else
			id = -1;
	}
	return id;
}

int Grid::traverseOverJoins3D(const int d1, const int d2, const int d3, const int dc1, const int dc2, const int dc3, int fromCellID) {
	int id = traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d1, dc1, fromCellID)));
	if (id == -1) {
		id = traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d3, dc3, fromCellID)));
		if (id == -1) {
			id = traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d2, dc2, id)));
			if (id == -1) {
				id = traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d3, dc3, fromCellID)));
				if (id == -1) {
					id = traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d1, dc1, fromCellID)));
					if (id == -1)
						id = traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d2, dc2, fromCellID)));
				}
			}
		}
	}
	return id;
}

int Grid::nextCell2D(const int plane, int fromCellID) {
	int id = fromCellID;
	Coords coords = unflatCoords(id);
	if ((int)(m_cells[id].xc[(plane+2)%3]+1)%2 != 0) {
		if (coords[(plane+1)%3] != coreCells[(plane+1)%3] - 1)
			id = traverse1D((plane+1)%3, 1, fromCellID);
		else
			id = -1;
	}
	else {
		if (coords[(plane+1)%3] != 0)
			id = traverse1D((plane+1)%3, -1, fromCellID);
		else
			id = -1;
	}
	if (id == -1 && coords[(plane+2)%3] != coreCells[(plane+2)%3] - 1)
		id = traverse1D((plane+2)%3, 1, fromCellID);
	return id;
}

int Grid::nextSnake(int fromCellID, int sourceCellID, const int dxc, const int dyc, const int dyz, int nd) {
	int id = traverse1D(0, dxc, fromCellID);
	if (id == -1 && nd > 1) {
		id = traverse1D(1, safe_round(m_cells[fromCellID].xc[1] - m_cells[sourceCellID].xc[1]) + dyc, sourceCellID);
		id = traverse1D(2, safe_round(m_cells[fromCellID].xc[2] - m_cells[sourceCellID].xc[2]), id);
	}
	if (id == -1 && nd > 2)
		id = traverse1D(2, safe_round(m_cells[fromCellID].xc[2] - m_cells[sourceCellID].xc[2]) + dyz, sourceCellID);
	return id;
}

int Grid::locate(int cx, int cy, int cz) {
	int cx_mod = cx - getLeftX();

	if (cx_mod >= getLeftX() && cx_mod < getRightX() && cy >= 0 && cy < coreCells[1] && cz >= 0 && cz < coreCells[2])
		return flatIndex(cx_mod, cy, cz);
	else
		return -1;
}

Coords Grid::nearestCoord(const Coords& original) {
	Coords nearest = original;
	nearest[0] = std::min(std::max(getLeftX(), nearest[0]), getRightX() - 1);
	nearest[1] = std::min(std::max(0, nearest[1]), coreCells[1] - 1);
	nearest[2] = std::min(std::max(0, nearest[2]), coreCells[2] - 1);
	return nearest;
}

bool Grid::withinGrid(const Coords& coords) {
	return coords[0] >= getLeftX() && coords[0] < getRightX() && coords[1] >= 0 && coords[1] < coreCells[1] && coords[2] >= 0 && coords[2] < coreCells[2];
}

bool Grid::joinExists(int dim, std::array<int, 3>& joinIDs) {
	return joinIDs[dim] >= 0 && joinIDs[dim] < (int)m_joins[dim].size();
}

Coords Grid::nextSnakeCoords(const Coords& fromCoords, const Coords& sourceCoords, const int dxc, const int dyc, const int dzc) {
	Coords nextCoords;
	if (fromCoords[0] + dxc >= getLeftX() && fromCoords[0] + dxc < getRightX()) {
		nextCoords[0] = fromCoords[0] + dxc;
		nextCoords[1] = fromCoords[1];
		nextCoords[2] = fromCoords[2];
	}
	else if (fromCoords[1] + dyc >= 0 && fromCoords[1] + dyc < coreCells[1]) {
		nextCoords[0] = sourceCoords[0];
		nextCoords[1] = fromCoords[1] + dyc;
		nextCoords[2] = fromCoords[2];
	}
	else if (fromCoords[2] + dzc >= 0 && fromCoords[2] + dzc < coreCells[2]) {
		nextCoords[0] = sourceCoords[0];
		nextCoords[1] = sourceCoords[1];
		nextCoords[2] = fromCoords[2] + dzc;
	}
	else {
		nextCoords[0] = -1;
		nextCoords[1] = -1;
		nextCoords[2] = -1;
	}

	return nextCoords;
}

Coords Grid::nextCausalCoords(const Coords& fromCoords, const Coords& sourceCoords) {
	int dir[3] = {0,0,0};
	for (int i = 0; i < 3; ++i){
		if (fromCoords[i] != sourceCoords[i])
			dir[i] = (fromCoords[i]-sourceCoords[i])/std::fabs(fromCoords[i]-sourceCoords[i]);
		else
			dir[i] = 1;
	}
	Coords beginCoords = sourceCoords;
	for (int i = 0; i < 3; ++i) {
		if (dir[i] < 0)
			beginCoords[i] -= 1;
	}

	Coords nextCoords = nextSnakeCoords(fromCoords, beginCoords, dir[0], dir[1], dir[2]);

	if (!withinGrid(nextCoords) && dir[0] == 1 && dir[1] == 1 && dir[2] == 1) {
		nextCoords = sourceCoords;
		nextCoords[0] -= 1;
		dir[0] = -1;
	}
	if (!withinGrid(nextCoords) && dir[0] == -1 && dir[1] == 1 && dir[2] == 1 && m_consts->nd > 1) {
		nextCoords = sourceCoords;
		nextCoords[1] -= 1;
		dir[0] = 1;
		dir[1] = -1;
	}
	if (!withinGrid(nextCoords) && dir[0] == 1 && dir[1] == -1 && dir[2] == 1 && m_consts->nd > 1) {
		nextCoords = sourceCoords;
		nextCoords[0] -= 1;
		nextCoords[1] -= 1;
		dir[0] = -1;
	}
	if (!withinGrid(nextCoords) && dir[0] == -1 && dir[1] == -1 && dir[2] == 1 && m_consts->nd > 2) {
		nextCoords = sourceCoords;
		nextCoords[2] -= 1;
		dir[0] = 1;
		dir[1] = 1;
		dir[2] = -1;
	}
	if (!withinGrid(nextCoords) && dir[0] == 1 && dir[1] == 1 && dir[2] == -1 && m_consts->nd > 2) {
		nextCoords = sourceCoords;
		nextCoords[0] -= 1;
		nextCoords[2] -= 1;
		dir[0] = -1;
	}
	if (!withinGrid(nextCoords) && dir[0] == -1 && dir[1] == 1 && dir[2] == -1 && m_consts->nd > 2) {
		nextCoords = sourceCoords;
		nextCoords[1] -= 1;
		nextCoords[2] -= 1;
		dir[0] = 1;
		dir[1] = -1;
	}
	if (!withinGrid(nextCoords) && dir[0] == 1 && dir[1] == -1 && dir[2] == -1 && m_consts->nd > 2) {
		nextCoords = sourceCoords;
		nextCoords[0] -= 1;
		nextCoords[1] -= 1;
		nextCoords[2] -= 1;
		dir[0] = -1;
	}

	return nextCoords;
}

int Grid::nextCausal(int fromCellID, int sourceCellID, int nd) {
	int dir[3] = {0,0,0};
	for (int i = 0; i < 3; ++i){
		if ((int)m_cells[fromCellID].xc[i] != (int)m_cells[sourceCellID].xc[i])
			dir[i] = (m_cells[fromCellID].xc[i]-m_cells[sourceCellID].xc[i])/std::fabs(m_cells[fromCellID].xc[i]-m_cells[sourceCellID].xc[i]);
		else
			dir[i] = 1;
	}
	int beginID = sourceCellID;
	if(dir[0] < 0 && dir[1] > 0 && dir[2] > 0)
		beginID = m_cells[sourceCellID].leftID[0];
	if(dir[0] > 0 && dir[1] < 0 && dir[2] > 0)
		beginID = m_cells[sourceCellID].leftID[1];
	if(dir[0] > 0 && dir[1] > 0 && dir[2] < 0)
		beginID = m_cells[sourceCellID].leftID[2];
	if(dir[0] < 0 && dir[1] < 0 && dir[2] > 0)
		beginID = m_cells[m_cells[sourceCellID].leftID[0]].leftID[1];
	if(dir[0] < 0 && dir[1] > 0 && dir[2] < 0)
		beginID = m_cells[m_cells[sourceCellID].leftID[0]].leftID[2];
	if(dir[0] > 0 && dir[1] < 0 && dir[2] < 0)
		beginID = m_cells[m_cells[sourceCellID].leftID[1]].leftID[2];
	if(dir[0] < 0 && dir[1] < 0 && dir[2] < 0)
		beginID = m_cells[m_cells[m_cells[sourceCellID].leftID[0]].leftID[1]].leftID[2];

	int id = nextSnake(fromCellID, beginID, dir[0], dir[1], dir[2], nd);

	if (id == -1 && dir[0] == 1 && dir[1] == 1 && dir[2] == 1) {
		id = traverse3D(0, 1, 2, -1, 0, 0, sourceCellID);
		dir[0] = -1;
	}
	if (id == -1 && dir[0] == -1 && dir[1] == 1 && dir[2] == 1 && nd > 1) {
		id = traverse3D(0, 1, 2, 0, -1, 0, sourceCellID);
		dir[0] = 1;
		dir[1] = -1;
	}
	if (id == -1 && dir[0] == 1 && dir[1] == -1 && dir[2] == 1 && nd > 1) {
		id = traverse3D(0, 1, 2, -1, -1, 0, sourceCellID);
		dir[0] = -1;
	}
	if (id == -1 && dir[0] == -1 && dir[1] == -1 && dir[2] == 1 && nd > 2) {
		id = traverse3D(0, 1, 2, 0, 0, -1, sourceCellID);
		dir[0] = 1;
		dir[1] = 1;
		dir[2] = -1;
	}
	if (id == -1 && dir[0] == 1 && dir[1] == 1 && dir[2] == -1 && nd > 2) {
		id = traverse3D(0, 1, 2, -1, 0, -1, sourceCellID);
		dir[0] = -1;
	}
	if (id == -1 && dir[0] == -1 && dir[1] == 1 && dir[2] == -1 && nd > 2) {
		id = traverse3D(0, 1, 2, 0, -1, -1, sourceCellID);
		dir[0] = 1;
		dir[1] = -1;
	}
	if (id == -1 && dir[0] == 1 && dir[1] == -1 && dir[2] == -1 && nd > 2) {
		id = traverse3D(0, 1, 2, -1, -1, -1, sourceCellID);
		dir[0] = -1;
	}
	return id;
}

void Grid::weakLink(const int dim, int lcellID, int rcellID) {
	if (lcellID != -1 && rcellID != -1) {
		m_joins[dim].emplace_back();
		m_cells[rcellID].ljoinID[dim] = m_joins[dim].size() - 1;
		m_cells[lcellID].rjoinID[dim] = m_joins[dim].size() - 1;

		GridJoin& join = m_joins[dim][m_joins[dim].size() - 1];

		join.rcellID = rcellID;
		join.lcellID = lcellID;

		for (int i = 0; i < 3; ++i)
			join.xj[i] = m_cells[join.rcellID].xc[i];
		join.xj[dim] = m_cells[join.rcellID].xc[dim] - 0.5;
		join.area = computeJoinArea(join.xj, dim, dx, geometry, m_consts->nd);
	}
}

void Grid::link(const int dim, int lcellID, int rcellID) {
	if (lcellID != -1 && rcellID != -1) {
		m_cells[rcellID].leftID[dim] = lcellID;
		m_cells[lcellID].rightID[dim] = rcellID;

		weakLink(dim, lcellID, rcellID);
	}
}

static int calcCoreCells(int nx, int nproc, int iproc) {
	int ncells = nx/nproc;
	int cells_left = nx - nproc*ncells;
	if (nproc - iproc <= cells_left)
		++ncells;
	return ncells;
}

static int calcLeftBoundaryPosition(int nx, int nproc, int iproc) {
	int start_xc = 0;
	for (int oproc = 0; oproc < iproc; ++oproc)
		start_xc += calcCoreCells(nx, nproc, oproc);
	return start_xc;
}

void Grid::initialise(std::shared_ptr<Constants> consts, const GridParameters& gp) {
	m_consts = std::move(consts);
	MPIW& mpihandler = MPIW::Instance();
	int rproc = mpihandler.getRank();
	int nproc = mpihandler.nProcessors();

	ncells = gp.ncells;
	sideLength = gp.sideLength;
	spatialOrder = gp.spatialOrder;
	geometry = m_consts->geometryParser.parseEnum(gp.geometry);
	for (int i = 0; i < m_consts->nd; ++i)
		dx[i] = gp.sideLength/(double)gp.ncells[0];

	coreCells = gp.ncells;
	coreCells[0] = calcCoreCells(gp.ncells[0], nproc, rproc);
	if (coreCells[0] == 0) {
		throw std::runtime_error("Grid::initialise: zero cells in processor (" +  std::to_string(mpihandler.getRank()) + ").");
	}

	// Set Grid wall locations.
	setLeftX(calcLeftBoundaryPosition(ncells[0], nproc, rproc));
	setRightX(getLeftX()+coreCells[0]);

	m_cellCollection.start("AllCells");
	m_cellCollection.start("GridCells");
	buildCells();
	m_cellCollection.stop("GridCells");

	for (GridCell& cell : m_cells)
		cell.vol = computeCellVolume(cell.xc[0], dx, geometry, m_consts->nd);

	// Parse boundary condition into enum.
	std::pair<std::array<Condition, 3>, std::array<Condition, 3>> leftRightBC;
	for (int i = 0; i < 3; ++i) {
		leftRightBC.first[i] = m_consts->conditionParser.parseEnum(gp.leftBC[i]);
		leftRightBC.second[i] = m_consts->conditionParser.parseEnum(gp.rightBC[i]);
	}

	// Initialise partition manager.
	if (mpihandler.nproc > 1)
		partition.initialise(coreCells[1]*coreCells[2]*(spatialOrder+1)*(UID::N + 1));

	// Build the boundaries.
	buildBoundaries(leftRightBC.first, leftRightBC.second);

	// Link the boundaries to the Grid.
	m_cellCollection.start("GhostCells");
	for (unsigned int dim = 0; dim < m_boundaries.size() / 2; ++dim) {
		m_cellCollection.start("GhostCells"+std::to_string(dim));
		boundaryLink(m_boundaries[2*dim + 0]);
		boundaryLink(m_boundaries[2*dim + 1]);
		m_cellCollection.stop("GhostCells"+std::to_string(dim));
	}
	m_cellCollection.stop("GhostCells");
	m_cellCollection.stop("AllCells");
	m_cellCollection.start("DeepGhostCells");
	MPIW::Instance().barrier();
	MPIW::Instance().serial([&] () {
		for (Bound& boundary : m_boundaries)
			boundaryLinkDeeper(boundary);
	});
	m_cellCollection.stop("DeepGhostCells");
}

void Grid::buildCells() {
	int proc_rank = MPIW::Instance().getRank();

	int ncells = coreCells[0] * coreCells[1] * coreCells[2];

	Logger::Instance().print<SeverityType::NOTICE>("Building Grid3D...\n");
	ProgressBar progBar(ncells, 1000);

	for (int index = 0; index < ncells; ++index){
		double dummy;

		m_cellCollection.add();
		Coords nextCoords = unflatCoords(index);
		for (int i = 0; i < 3; ++i)
			m_cells[index].xc[i] = nextCoords[i];
		for (int i = 0; i < m_consts->nd; ++i)
			m_cells[index].xc[i] += 0.5;
		m_cells[index].xc[0] += getLeftX();

		progBar.update(index);
		if (progBar.timeToUpdate()) {
			progBar.update(index);
			Logger::Instance().print<SeverityType::INFO>(progBar.getFullString(), "\r");
		}
	}

	for (int i = 0; i < coreCells[0]; ++i) {
		for (int j = 0; j < coreCells[1]; ++j) {
			for (int k = 0; k < coreCells[2]; ++k) {
				if (i > 0)
					link(0, flatIndex(i-1, j, k), flatIndex(i, j, k));
				if (j > 0)
					link(1, flatIndex(i, j-1, k), flatIndex(i, j, k));
				if (k > 0)
					link(2, flatIndex(i, j, k-1), flatIndex(i, j, k));
			}
		}
	}

	progBar.end();
	Logger::Instance().print<SeverityType::NOTICE>(progBar.getFinalString(), '\n');
}

void Grid::buildCausal(const Coords& sourceCoords) {
	Coords startCoords = nearestCoord(sourceCoords);

	for (Coords c = startCoords; withinGrid(c); c = nextCausalCoords(c, startCoords)) {
		m_causalIndices.push_back(flatIndex(c[0] - getLeftX(), c[1], c[2]));
	}
}

void Grid::buildBoundaries(const std::array<Condition, 3>& leftBC, const std::array<Condition, 3>& rightBC) {
	MPIW& mpihandler = MPIW::Instance();
	int nproc = mpihandler.nProcessors();
	int rproc = mpihandler.getRank();

	int left_rank = (rproc - 1 < 0) ? rproc - 1 + nproc : rproc - 1;
	int right_rank = (rproc + 1 >= nproc) ? rproc + 1 - nproc : rproc + 1;
	bool periodic_and_mpi = nproc > 1 && leftBC[0] == Condition::PERIODIC && rightBC[0] == Condition::PERIODIC;

	if (m_consts->nd > 0) {
		if (mpihandler.getRank() != 0 || periodic_and_mpi)
			m_boundaries.push_back(Bound(0, Condition::PARTITION, left_rank));
		else
			m_boundaries.push_back(Bound(0, leftBC[0]));
		if (rproc != nproc - 1 || periodic_and_mpi)
			m_boundaries.push_back(Bound(3, Condition::PARTITION, right_rank));
		else
			m_boundaries.push_back(Bound(3, rightBC[0]));
	}
	if (m_consts->nd > 1) {
		m_boundaries.push_back(Bound(1, leftBC[1]));
		m_boundaries.push_back(Bound(4, rightBC[1]));
	}
	if (m_consts->nd > 2) {
		m_boundaries.push_back(Bound(2, leftBC[2]));
		m_boundaries.push_back(Bound(5, rightBC[2]));
	}
}

void Grid::boundaryLink(Bound& boundary) {
	int dim = boundary.face%3;
	bool isLeft = boundary.face < 3;
	int linkCellID = isLeft ? 0 : flatIndex(dim == 0 ? coreCells[0]-1 : 0, dim == 1 ? coreCells[1]-1 : 0, dim == 2 ? coreCells[2]-1 : 0);

	if (boundary.condition == Condition::PARTITION)
		m_cellCollection.start(isLeft ? "LeftPartitionCells" : "RightPartitionCells");
	while (linkCellID != -1) {
		int ghostCellID = m_cellCollection.add();
		GridCell& ghostCell = m_cells[ghostCellID];
		boundary.ghostCellIDs.push_back(ghostCellID);

		for (int idim = 0; idim < 3; ++idim)
			ghostCell.xc[idim] = m_cells[linkCellID].xc[idim];
		ghostCell.xc[dim] += isLeft ? -1 : 1;

		link(dim, isLeft ? ghostCellID : linkCellID, isLeft ?  linkCellID : ghostCellID);

		linkCellID = nextCell2D(dim, linkCellID);
	}
	if (boundary.condition == Condition::PARTITION)
		m_cellCollection.stop(isLeft ? "LeftPartitionCells" : "RightPartitionCells");
}

void Grid::boundaryLinkDeeper(Bound& boundary ) {
	int dim = boundary.face%3;
	bool isLeft = boundary.face < 3;

	for (int ghostCellID : boundary.ghostCellIDs) {
		int currGhostID = ghostCellID;

		for (int ig = 0; ig < spatialOrder; ++ig) {
			int oldGhostID = currGhostID;
			currGhostID = m_cellCollection.add();
			m_cells[currGhostID].xc[0] = -100;

			if (isLeft)
				m_cells[oldGhostID].leftID[dim] = currGhostID;
			else
				m_cells[oldGhostID].rightID[dim] = currGhostID;
		}
	}
}

void Grid::applyBCs() {
	for (Bound& boundary : m_boundaries) {
		int dim = boundary.face%3;

		switch(boundary.condition) {
			case(Condition::FREE):
				if (boundary.face < 3) {
					for (int ghostCellID : boundary.ghostCellIDs) {
						GridCell& cell = m_cells[right(dim, ghostCellID)];

						for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = left(dim, currGhostID)) {
							GridCell& ghost = m_cells[currGhostID];
							for(int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
						}
					}
				}
				else {
					for (int ghostCellID : boundary.ghostCellIDs) {
						GridCell& cell = m_cells[left(dim, ghostCellID)];

						for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = right(dim, currGhostID)) {
							GridCell& ghost = m_cells[currGhostID];
							for(int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
						}
					}
				}
				break;
			case(Condition::INFLOW):
				if (boundary.face < 3) {
					for (int ghostCellID : boundary.ghostCellIDs) {
						GridCell& cell = m_cells[right(dim, ghostCellID)];

						for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = left(dim, currGhostID)) {
							GridCell& ghost = m_cells[currGhostID];
							for(int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];
							if (cell.Q[UID::VEL+dim] < 0)
								ghost.Q[UID::VEL+dim] = -1.0*cell.Q[UID::VEL+dim];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
						}
					}
				}
				else {
					for (int ghostCellID : boundary.ghostCellIDs) {
						GridCell& cell = m_cells[left(dim, ghostCellID)];

						for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = right(dim, currGhostID)) {
							GridCell& ghost = m_cells[currGhostID];
							for (int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];
							if (cell.Q[UID::VEL+dim] > 0)
								ghost.Q[UID::VEL+dim] = -1.0*cell.Q[UID::VEL+dim];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
						}
					}
				}
				break;
			case(Condition::OUTFLOW):
				if (boundary.face < 3) {
					for (int ghostCellID : boundary.ghostCellIDs) {
						GridCell& cell = m_cells[right(dim, ghostCellID)];

						for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = left(dim, currGhostID)) {
							GridCell& ghost = m_cells[currGhostID];
							for (int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];
							if (cell.Q[UID::VEL+dim] > 0)
								ghost.Q[UID::VEL+dim] = -1.0*cell.Q[UID::VEL+dim];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
						}
					}
				}
				else {
					for (int ghostCellID : boundary.ghostCellIDs) {
						GridCell& cell = m_cells[left(dim, ghostCellID)];

						for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = right(dim, currGhostID)) {
							GridCell& ghost = m_cells[currGhostID];
							for (int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];
							if(cell.Q[UID::VEL+dim] < 0)
								ghost.Q[UID::VEL+dim] = -1.0*cell.Q[UID::VEL+dim];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
						}
					}
				}
				break;
			case(Condition::PERIODIC):
				if (boundary.face < 3) {
					for (int ghostCellID : boundary.ghostCellIDs) {
						Coords cellCoords = unflatCoords(right(dim, ghostCellID));
						cellCoords[dim] += coreCells[dim] - 1;

						int currCellID = flatIndex(cellCoords[0], cellCoords[1], cellCoords[2]);
						for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = left(dim, currGhostID)) {
							GridCell& ghost = m_cells[currGhostID];
							GridCell& cell = m_cells[currCellID];
							for(int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
							currCellID = left(dim, currCellID);
						}
					}
				}
				else {
					for (int ghostCellID : boundary.ghostCellIDs) {
						Coords cellCoords = unflatCoords(left(dim, ghostCellID));
						cellCoords[dim] -= coreCells[dim] - 1;

						int currCellID = flatIndex(cellCoords[0], cellCoords[1], cellCoords[2]);
						for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = right(dim, currGhostID)) {
							GridCell& ghost = m_cells[currGhostID];
							GridCell& cell = m_cells[currCellID];
							for(int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
							currCellID = right(dim, currCellID);
						}
					}
				}
				break;
			case(Condition::REFLECTING):
				if (boundary.face < 3) {
					for (int ghostCellID : boundary.ghostCellIDs) {
						for (int currGhostID = ghostCellID, currCellID = right(dim, ghostCellID); currGhostID != -1; currGhostID = left(dim, currGhostID)) {
							currCellID = right(dim, currCellID);
							GridCell& ghost = m_cells[currGhostID];
							GridCell& cell = m_cells[currCellID];

							for (int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];

							ghost.Q[UID::VEL+dim] = -cell.Q[UID::VEL+dim];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
						}
					}
				}
				else {
					for (int ghostCellID : boundary.ghostCellIDs) {
						for (int currGhostID = ghostCellID, currCellID = left(dim, ghostCellID); currGhostID != -1; currGhostID = right(dim, currGhostID)) {
							currCellID = left(dim, currCellID);
							GridCell& ghost = m_cells[currGhostID];
							GridCell& cell = m_cells[currCellID];
							for(int iu = 0; iu < UID::N; iu++)
								ghost.Q[iu] = cell.Q[iu];
							ghost.Q[UID::VEL+dim] = -cell.Q[UID::VEL+dim];
							ghost.heatCapacityRatio = cell.heatCapacityRatio;
						}
					}
				}
				break;
			case(Condition::PARTITION):
				MPIW& mpihandler = MPIW::Instance();

				partition.resetBuffer();

				for (int ghostCellID : boundary.ghostCellIDs) {
					int currCellID = ghostCellID;
					for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = boundary.face < 3 ? left(dim, currGhostID) : right(dim, currGhostID)) {
						currCellID = boundary.face < 3 ? right(dim, currCellID) : left(dim, currCellID);
						GridCell& cell = m_cells[currCellID];
						for (int iu = 0; iu < UID::N; ++iu)
							partition.addSendItem(cell.Q[iu]);
						partition.addSendItem(cell.heatCapacityRatio);
					}
				}

				bool isPeriodic = (mpihandler.getRank() == 0 && boundary.face == 0) || (mpihandler.getRank() == mpihandler.nProcessors()-1 && boundary.face == 3);
				SendID tag = isPeriodic ? SendID::PERIODIC_MSG : SendID::PARTITION_MSG;

				partition.exchangeData(boundary.targetProcessor, tag);

				for (int ghostCellID : boundary.ghostCellIDs) {
					for (int currGhostID = ghostCellID; currGhostID != -1; currGhostID = (boundary.face < 3) ? left(dim, currGhostID) : right(dim, currGhostID)) {
						GridCell& ghost = m_cells[currGhostID];
						for (int iu = 0; iu < UID::N; ++iu)
							ghost.Q[iu] = partition.getRecvItem();
						ghost.heatCapacityRatio = partition.getRecvItem();
					}
				}
				break;
		}
	}
}

/**
 * @brief Calculates the area between two GridCell objects given its location and the geometry of Grid3D.
 * @param xj
 * Coordinates of the GridJoin between the two GridCell objects.
 * @param dim
 * Direction normal to the face of the GridJoin.
 * @return The area.
 */
double Grid::computeJoinArea(const Vec3& xj, const int dim, const Vec3& dx, Geometry geometry, int nd) {
	double PI = 3.14159265359;
	double area = 0, r1, r2;
	if (geometry == Geometry::CARTESIAN) {
		area = 1.0;
		for (int i = 0; i < nd; ++i) {
			if (i != dim)
				area *= dx[i];
		}
	}
	else if (geometry == Geometry::CYLINDRICAL){
		/* Cylindrical polars (ND <= 2) */
		if (dim == 0) {
			r1 = xj[0]*dx[0];
			area = 2.0*PI*r1;
			if (nd >= 2)
				area *= dx[1];
		}
		else if (dim == 1) {
			r1 = xj[0]*dx[0] - 0.5*dx[0];
			r2 = xj[0]*dx[0] + 0.5*dx[0];
			area = PI*(r2 - r1)*(r2 + r1);
		}
	}
	else if (geometry == Geometry::SPHERICAL) {
		/* Spherical polars (ND = 1 only) */
		r1 = xj[0]*dx[0];
		area = 4.0*PI*r1*r1;
	}
	return area;
}

/**
 * @brief Calculates the volume of a GridCell given its location and the geometry of Grid3D.
 * @param cptr
 * The GridCell object for which the volume will be calculated.
 */
double Grid::computeCellVolume(double rc, const Vec3& dx, Geometry geometry, int nd) {
	double PI = 3.14159265359;
	double volume = 0, r1, r2;
	if (geometry == Geometry::CARTESIAN) {
		/* Cartesian */
		volume = 1.0;
		for (int i = 0; i < nd; ++i)
			volume *= dx[i];
	}
	else if (geometry == Geometry::CYLINDRICAL) {
		/* Cylindrical polars (ND <= 2) */
		r2 = rc*dx[0] + 0.5*dx[0];
		r1 = rc*dx[0] - 0.5*dx[0];
		volume = PI*(r2 - r1)*(r2 + r1);
		if (nd >= 2)
			volume *= dx[1];
	}
	else if (geometry == Geometry::SPHERICAL) {
		/* Spherical polars (ND = 1 only) */
		r2 = rc*dx[0] + 0.5*dx[0];
		r1 = rc*dx[0] - 0.5*dx[0];
		volume = 4.0*PI*(r2 - r1)*(r2*r2 + r1*r2 + r1*r1)/3.0;
	}
	return volume;
}

int Grid::getRayPlane(const Vec3& xc, const Vec3& xs) const {
	int plane = 0;
	double dxcs_max = 0;
	for (int i = 0; i < m_consts->nd; ++i) {
		double dxcs = fabs(xc[i] - xs[i]);
		if (dxcs > dxcs_max) {
			dxcs_max = dxcs;
			plane = i;
		}
	}
	return plane;
}

void Grid::calculateNearestNeighbours(const std::array<double, 3>& star_pos) {
	for (int index = 0; index < coreCells[0]*coreCells[1]*coreCells[2]; ++index) {

		int plane = getRayPlane(m_cells[index].xc, star_pos);
		int irot[3] = {(plane+1)%3, (plane+2)%3, (plane%3)};
		double d[3] = {0.0, 0.0, 0.0};
		for(int i = 0; i < 3; i++)
			d[i] = m_cells[index].xc[irot[i]] - star_pos[irot[i]];
		int s[3] = {d[0] < -1.0/10.0 ? -1 : 1, d[1] < -1.0/10.0 ? -1 : 1, d[2] < -1.0/10.0 ? -1 : 1};
		int LR[3] = {std::abs(d[0]) < 1.0/10.0 ? 0 : s[0], std::abs(d[1]) < 1.0/10.0 ? 0 : s[1], std::abs(d[2]) < 1.0/10.0 ? 0 : s[2]};
		m_cells[index].neighbourIDs[0] = traverse3D(irot[0], irot[1], irot[2], 0, 0, -LR[2], index);
		m_cells[index].neighbourIDs[1] = traverse3D(irot[0], irot[1], irot[2], 0, -LR[1], -LR[2], index);
		m_cells[index].neighbourIDs[2] = traverse3D(irot[0], irot[1], irot[2], -LR[0], 0, -LR[2], index);
		m_cells[index].neighbourIDs[3] = traverse3D(irot[0], irot[1], irot[2], -LR[0], -LR[1], -LR[2], index);
		if (m_cells[index].neighbourIDs[0] == -1)
			m_cells[index].neighbourIDs[0] = Grid::traverseOverJoins3D(irot[0], irot[1], irot[2], 0, 0, -LR[2], index);
		if (m_cells[index].neighbourIDs[1] == -1)
			m_cells[index].neighbourIDs[1] = Grid::traverseOverJoins3D(irot[0], irot[1], irot[2], 0, -LR[1], -LR[2], index);
		if (m_cells[index].neighbourIDs[2] == -1)
			m_cells[index].neighbourIDs[2] = Grid::traverseOverJoins3D(irot[0], irot[1], irot[2], -LR[0], 0, -LR[2], index);
		if (m_cells[index].neighbourIDs[3] == -1)
			m_cells[index].neighbourIDs[3] = Grid::traverseOverJoins3D(irot[0], irot[1], irot[2], -LR[0], -LR[1], -LR[2], index);

		double ic[3] = {(int)m_cells[index].xc[irot[0]]-0.5*(s[2]*d[0]/d[2]),	(int)m_cells[index].xc[irot[1]]-0.5*(s[2]*d[1]/d[2]),	(int)m_cells[index].xc[irot[2]]-0.5*(s[2])};
		double delta[2] = {std::abs(2.0*ic[0]-2.0*(int)m_cells[index].xc[irot[0]]+s[0]), std::abs(2.0*ic[1]-2.0*(int)m_cells[index].xc[irot[1]]+s[1])};
		m_cells[index].neighbourWeights[0] = (std::abs(d[2]) > 0.9) ? delta[0]*delta[1] : 0;
		m_cells[index].neighbourWeights[1] = ((std::abs(d[1]) > 0.9) && (std::abs(d[2]) > 0.9)) ? delta[0]*(1.0-delta[1]) : 0;
		m_cells[index].neighbourWeights[2] = ((std::abs(d[0]) > 0.9) && (std::abs(d[2]) > 0.9)) ? (1.0-delta[0])*delta[1] : 0;
		m_cells[index].neighbourWeights[3] = ((std::abs(d[0]) > 0.9) && (std::abs(d[1]) > 0.9) && (std::abs(d[2]) > 0.9)) ? (1.0-delta[0])*(1.0-delta[1]) : 0;
	}
}

Bound::Bound(int face, const Condition bcond, int target_proc)
: face(face)
, condition(bcond)
, targetProcessor(target_proc)
{
}
