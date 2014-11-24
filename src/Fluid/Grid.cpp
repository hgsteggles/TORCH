#include "Grid.hpp"
#include "Boundary.hpp"

#include <cmath>


static int safe_round(double val) {
	if (val < 0)
		return (int)(val-0.5);
	else
		return (int)(val+0.5);
}

Grid::Grid()
{ }

Grid::Grid(Grid&& o)
{
	ncells = o.ncells;
	coreCells = o.coreCells;
	m_leftX = o.m_leftX;
	m_rightX = o.m_rightX;
	currentTime = o.currentTime;
	deltatime = o.deltatime;
	dx = o.dx;
	geometry = o.geometry;
	spatialOrder = o.spatialOrder;
	sideLength = o.sideLength;
	std::swap(m_cells, o.m_cells);
	for (int i = 0; i < 3; ++i) {
		std::swap(m_joins[i], o.m_joins[i]);
		std::swap(m_leftBoundaries[i], o.m_leftBoundaries[i]);
		std::swap(m_rightBoundaries[i], o.m_rightBoundaries[i]);
	}
}

Grid& Grid::operator=(Grid&& o) {
	ncells = o.ncells;
	coreCells = o.coreCells;
	m_leftX = o.m_leftX;
	m_rightX = o.m_rightX;
	currentTime = o.currentTime;
	deltatime = o.deltatime;
	dx = o.dx;
	geometry = o.geometry;
	spatialOrder = o.spatialOrder;
	sideLength = o.sideLength;
	std::swap(m_cells, o.m_cells);
	for (int i = 0; i < 3; ++i) {
		std::swap(m_joins[i], o.m_joins[i]);
		std::swap(m_leftBoundaries[i], o.m_leftBoundaries[i]);
		std::swap(m_rightBoundaries[i], o.m_rightBoundaries[i]);
	}
	return *this;
}

Grid::~Grid() {
	m_cells.getAll().dispose();

	for (JoinContainer& join_dim : m_joins) {
		join_dim.dispose();
	}
	for (Boundary* bptr : m_leftBoundaries) {
		if (bptr != nullptr) {
			delete bptr;
			bptr = nullptr;
		}
	}
	for (Boundary* bptr : m_rightBoundaries) {
		if (bptr != nullptr) {
			delete bptr;
			bptr = nullptr;
		}
	}
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

CellContainer& Grid::getCells() {
	return m_cells.getAll();
}

CellContainer& Grid::getWindCells() {
	return m_cells.getFirst();
}

CellContainer& Grid::getCausalCells() {
	return m_cells.getSecond();
}

JoinContainer& Grid::getJoins(int dim) {
	return m_joins[dim];
}

std::array<Boundary*, 3>& Grid::getLeftBoundaries() {
	return m_leftBoundaries;
}

std::array<Boundary*, 3>& Grid::getRightBoundaries() {
	return m_rightBoundaries;
}

const CellContainer& Grid::getCells() const {
	return m_cells.getAll();
}

const CellContainer& Grid::getWindCells() const {
	return m_cells.getFirst();
}

const CellContainer& Grid::getCausalCells() const {
	return m_cells.getSecond();
}

const JoinContainer& Grid::getJoins(int dim) const {
	return m_joins[dim];
}

const std::array<Boundary*, 3>& Grid::getLeftBoundaries() const {
	return m_leftBoundaries;
}

const std::array<Boundary*, 3>& Grid::getRightBoundaries() const {
	return m_rightBoundaries;
}

/**
 * @brief Traverses over GridCells via GridCell::right and GridCell::left in 1D.
 * @param dim
 * The dimension to traverse along (0, 1, 2).
 * @param dxc
 * No. of GridCells to traverse over. Passing 1 would return the adjacent GridCell.
 * @param cptr
 * GridCell to start traversing from.
 * @return The GridCell traversed to. Will return <code> nullptr </code> if it doesn't exist.
 */
GridCell* Grid::traverse1D(const int dim, const int dxc, GridCell* cptr) {
	GridCell* newcptr = cptr;
	int adxc = std::abs(dxc);
	for (int i = 0; i < adxc && newcptr != nullptr; ++i)
		newcptr = (dxc > 0) ? newcptr->right[dim] : newcptr->left[dim];
	return newcptr;
}

/**
 * @brief Traverses over GridCells via GridCell::right and GridCell::left in 3D.
 * blach
 * @param d1
 * Direction 1.
 * @param d2
 * Direction 2.
 * @param d3
 * Direction 3.
 * @param dc1
 * No. of GridCells to traverse along direction d1.
 * @param dc2
 * No. of GridCells to traverse along direction d2.
 * @param dc3
 * No. of GridCells to traverse along direction d3.
 * @param cptr
 * Pointer to the GridCell to start from.
 * @return The GridCell traversed to. Will return <code> nullptr </code> if it doesn't exist.
 */
GridCell* Grid::traverse3D(const int d1, const int d2, const int d3, const int dc1, const int dc2, const int dc3, GridCell* cptr) {
	return traverse1D(d3, dc3, traverse1D(d2, dc2, traverse1D(d1, dc1, cptr)));
}

/**
 * @brief Traverses over GridCells via GridCell::rjoin and GridCell::ljoin in 1D.
 * @param dim
 * The dimension to traverse along (0, 1, 2).
 * @param dc
 * No. of GridCells to traverse over. Passing 1 would return the adjacent GridCell.
 * @param cptr
 * Pointer to the GridCell to start from.
 * @return The GridCell traversed to. Will return <code> nullptr </code> if it doesn't exist.
 */
GridCell* Grid::traverseOverJoins1D(const int dim, const int dxc, GridCell* cptr) {
	GridCell* newcptr = cptr;
	for(int i = 0; i < abs(dxc) && newcptr != nullptr; ++i){
		if (dxc > 0 && newcptr->rjoin[dim] != nullptr)
			newcptr = newcptr->rjoin[dim]->rcell;
		else if (dxc < 0 && newcptr->ljoin[dim] != nullptr)
			newcptr = newcptr->ljoin[dim]->lcell;
		else
			newcptr = nullptr;
	}
	return newcptr;
}

/**
 * @brief Traverses over GridCells via GridCell::rjoin and GridCell::ljoin in 1D.
 * @param d1
 * Direction 1.
 * @param d2
 * Direction 2.
 * @param d3
 * Direction 3.
 * @param dc1
 * No. of GridCells to traverse along direction d1.
 * @param dc2
 * No. of GridCells to traverse along direction d2.
 * @param dc3
 * No. of GridCells to traverse along direction d3.
 * @param cptr
 * Pointer to the GridCell to start traversing from.
 * @return The GridCell traversed to. Will return <code> nullptr </code> if it doesn't exist.
 */
GridCell* Grid::traverseOverJoins3D(const int d1, const int d2, const int d3, const int dc1, const int dc2, const int dc3, GridCell* cptr) {
	GridCell* newcptr = traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d1, dc1, cptr)));
	if (newcptr == nullptr)
		newcptr = traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d3, dc3, cptr)));
	if (newcptr == nullptr)
		newcptr = traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d2, dc2, cptr)));
	if (newcptr == nullptr)
		newcptr = traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d3, dc3, cptr)));
	if (newcptr == nullptr)
		newcptr = traverseOverJoins1D(d2, dc2, traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d1, dc1, cptr)));
	if (newcptr == nullptr)
		newcptr = traverseOverJoins1D(d3, dc3, traverseOverJoins1D(d1, dc1, traverseOverJoins1D(d2, dc2, cptr)));
	return newcptr;
}

/**
 * @brief Provides the next GridCell along a path that traverses a 2D plane which includes the passed GridCell.
 * @param plane
 * The dimension normal to the plane.
 * @param cptr
 * Pointer to previous GridCell.
 * @return Pointer to next GridCell. If there are no more GridCells on this path then <code> nullptr </code> is returned.
 */
GridCell* Grid::nextCell2D(const int plane, GridCell* cptr) {
	GridCell* newcptr = nullptr;
	if ((int)(cptr->xc[(plane+2)%3]+1)%2 != 0)
		newcptr = traverse1D((plane+1)%3, 1, cptr);
	else
		newcptr = traverse1D((plane+1)%3, -1, cptr);
	if (newcptr == nullptr)
		newcptr = traverse1D((plane+2)%3, 1, cptr);
	return newcptr;
}

/**
 * @brief Provides next GridCell along a specific path.
 * The path starts from srcptr and traverses along the positive or negative x direction (sign of dxc) until a
 * <code> nullptr </code> pointer is encountered. The next GridCell in this case would be along the y direction
 * (specified by dyc) from a GridCell that has the x and z coordinates of the srcptr and the y coordinate of the
 * previous GridCell. At some point traversing in the y direction will encounter a nullptr pointer a so step is
 * taken in the z direction starting from a GridCell with x and y coordinates of the GridCell pointed to
 * by srcptr and z coordinate of the previous GridCell. E.g. (dxc,dyc,dzc)=(1,1,1) leads to all GridCell with
 * x,y,z coordinate equal to or greater than the GridCell pointed to by srcptr being sampled along this route.
 *
 * @param srcptr
 * The GridCell that traversal in all directions starts from.
 * @param cptr
 * Pointer to previous GridCell.
 * @param dxc
 * GridCell to traverse in the x direction at every step.
 * @param dyc
 * GridCells to traverse in the y direction at every step.
 * @param dyz
 * GridCells to traverse in the z direction at every step.
 * @return Pointer to next GridCell. If there are no more GridCells on this path then <code> nullptr </code> is returned.
 */
GridCell* Grid::nextSnake(GridCell* cptr, GridCell* srcptr, const int dxc, const int dyc, const int dyz, int nd) {
	GridCell* newcptr = traverse1D(0, dxc, cptr);
	if (newcptr == nullptr && nd > 1) {
		newcptr = traverse1D(1, safe_round(cptr->xc[1] - srcptr->xc[1]) + dyc, srcptr);
		newcptr = traverse1D(2, safe_round(cptr->xc[2] - srcptr->xc[2]), newcptr);
	}
	if (newcptr == nullptr && nd > 2)
		newcptr = traverse1D(2, safe_round(cptr->xc[2] - srcptr->xc[2]) + dyz, srcptr);
	return newcptr;
}

/**
 * @brief Provides next GridCell along a causal path from a specified GridCell.
 * @param cptr
 * Pointer to previous GridCell.
 * @param srcptr
 * Pointer to the first GridCell along the path.
 * @return Pointer to next GridCell. If there are no more GridCells on this path then <code> nullptr </code> is returned.
 */
GridCell* Grid::nextCausal(GridCell* cptr, GridCell* srcptr, int nd) {
	int dir[3] = {0,0,0};
	for (int i = 0; i < 3; ++i){
		if ((int)cptr->xc[i] != (int)srcptr->xc[i])
			dir[i] = (cptr->xc[i]-srcptr->xc[i])/std::fabs(cptr->xc[i]-srcptr->xc[i]);
		else
			dir[i] = 1;
	}
	GridCell* beginptr = srcptr;
	if(dir[0] < 0 && dir[1] > 0 && dir[2] > 0)
		beginptr = srcptr->left[0];
	if(dir[0] > 0 && dir[1] < 0 && dir[2] > 0)
		beginptr = srcptr->left[1];
	if(dir[0] > 0 && dir[1] > 0 && dir[2] < 0)
		beginptr = srcptr->left[2];
	if(dir[0] < 0 && dir[1] < 0 && dir[2] > 0)
		beginptr = srcptr->left[0]->left[1];
	if(dir[0] < 0 && dir[1] > 0 && dir[2] < 0)
		beginptr = srcptr->left[0]->left[2];
	if(dir[0] > 0 && dir[1] < 0 && dir[2] < 0)
		beginptr = srcptr->left[1]->left[2];
	if(dir[0] < 0 && dir[1] < 0 && dir[2] < 0)
		beginptr = srcptr->left[0]->left[1]->left[2];

	GridCell* newcptr = nextSnake(cptr, beginptr, dir[0], dir[1], dir[2], nd);

	if(newcptr == nullptr && dir[0] == 1 && dir[1] == 1 && dir[2] == 1) {
		newcptr = traverse3D(0, 1, 2, -1, 0, 0, srcptr);
		dir[0] = -1;
	}
	if(newcptr == nullptr && dir[0] == -1 && dir[1] == 1 && dir[2] == 1 && nd > 1) {
		newcptr = traverse3D(0, 1, 2, 0, -1, 0, srcptr);
		dir[0] = 1;
		dir[1] = -1;
	}
	if(newcptr == nullptr && dir[0] == 1 && dir[1] == -1 && dir[2] == 1 && nd > 1) {
		newcptr = traverse3D(0, 1, 2, -1, -1, 0, srcptr);
		dir[0] = -1;
	}
	if(newcptr == nullptr && dir[0] == -1 && dir[1] == -1 && dir[2] == 1 && nd > 2) {
		newcptr = traverse3D(0, 1, 2, 0, 0, -1, srcptr);
		dir[0] = 1;
		dir[1] = 1;
		dir[2] = -1;
	}
	if(newcptr == nullptr && dir[0] == 1 && dir[1] == 1 && dir[2] == -1 && nd > 2) {
		newcptr = traverse3D(0, 1, 2, -1, 0, -1, srcptr);
		dir[0] = -1;
	}
	if(newcptr == nullptr && dir[0] == -1 && dir[1] == 1 && dir[2] == -1 && nd > 2) {
		newcptr = traverse3D(0, 1, 2, 0, -1, -1, srcptr);
		dir[0] = 1;
		dir[1] = -1;
	}
	if(newcptr == nullptr && dir[0] == 1 && dir[1] == -1 && dir[2] == -1 && nd > 2) {
		newcptr = traverse3D(0, 1, 2, -1, -1, -1, srcptr);
		dir[0] = -1;
	}
	return newcptr;
}

/**
 * @brief Locates GridCell at the specified grid coordinates.
 * @param x
 * x grid coordinate.
 * @param y
 * y grid coordinate.
 * @param z
 * z grid coordinate.
 * @return The located GridCell if it exists. If not then <code> nullptr </code>.
 */
GridCell* Grid::locate(const Coords& x, GridCell* const fromCell) {
	GridCell* cptr = fromCell;
	for (int dim = 0; dim < 3; ++dim) {
		for (int i = 0; (cptr != nullptr) && ((int)(cptr->xc[dim]) != x[dim]); ++i)
			cptr = cptr->right[dim];
	}
	return cptr;
}

/**
 * @brief Locates the GridCell nearest to the passed coordinates given a Grid's first and last GridCell.
 * @param x
 * @param firstLastCell
 * @return The located GridCell if it exists. If not then <code> nullptr </code>.
 */
GridCell* Grid::getNearestCell(const Coords& x, const GridCellPair& firstLastCell) {
	Coords search_x = x;
	if (x[0] < safe_round(firstLastCell.first->xc[0]))
		search_x[0] = firstLastCell.first->xc[0];
	else if (x[0] > safe_round(firstLastCell.second->xc[0]))
		search_x[0] = firstLastCell.second->xc[0];
	if (x[1] < safe_round(firstLastCell.first->xc[1]))
		search_x[1] = firstLastCell.first->xc[1];
	else if (x[1] > safe_round(firstLastCell.second->xc[1]))
		search_x[1] = firstLastCell.second->xc[1];
	if (x[2] < safe_round(firstLastCell.first->xc[2]))
		search_x[2] = firstLastCell.first->xc[2];
	else if (x[2] > safe_round(firstLastCell.second->xc[2]))
		search_x[2] = firstLastCell.second->xc[2];
	return locate(search_x, firstLastCell.first);
}

