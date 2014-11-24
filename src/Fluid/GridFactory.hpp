/**
 * Provides the GridFactory class.
 * @file GridFactory.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef GRIDBUILDER_HPP_
#define GRIDBUILDER_HPP_

#include "Container.hpp"
#include "Common.hpp"
#include "Grid.hpp"

#include <array>
#include <vector>

class GridCell;
class GridJoin;
class Boundary;
class GridParameters;

/**
 * @class GridFactory
 *
 * @brief The GridFactory builds a doubly linked 3D grid and the Boundary's that enclose it.
 *
 * A GridFactory generates GridCells and links them to their nearest neighbours via direct pointers and via GridJoins which each have
 * pointers to the GridCells either side of them.
 * It then builds Boundary's that link, only via GridJoins (for automatic Boundary detection), to the faces of the Grid.
 *
 * GridCell* <==> GridJoin* <==> GridCell* <==> GridJoin <==> GridCell*(ghost)      GridCell*(ghost)...
 *           <=================>                                               <==>
 *
 * Finally, the Grid is built by adding GridCells to a DualIntrusiveContainer (in causal order). A DualIntrusiveContainer is needed
 * because the Grid is comprised of GridCells in a wind region and GridCells outside a wind region which may need to be dealt with
 * separately. DualIntrusiveContainer allows three types of iteration: over wind GridCells, over other GridCells and over all GridCells.
 *
 * During the creation, GridCell volumes and GridJoin areas are calculated.
 *
 * @see Grid
 * @see GridCell
 * @see GridJoin
 * @see ExternalBoundary
 * @see Partition
 * @see IntrusiveContainer
 * @see DualIntrusiveContainer
 *
 * @version 0.8, 24/11/2014
 */
class GridFactory {
	using CellContainerPair = std::pair<CellContainer, CellContainer>;
	using BoundaryVecPair = std::pair<BoundaryVec, BoundaryVec>;
	using ConditionVecPair = std::pair<ConditionVec, ConditionVec>;
	using GridCellPair = std::pair<GridCell* const, GridCell* const>;
public:
	static Grid createGrid(const std::shared_ptr<Constants>& c, const GridParameters& gp, const std::array<int, 3>& star_pos, int radius);
private:

	static GridCellPair buildGrid(const Coords& ncells, int nd, int start_xc, std::array<JoinContainer, 3>& joins);
	static DualIntrusiveContainer<GridCell> zipList(const GridCellPair& firstLastCell, const Coords& xs, int radius, int nd);
	static BoundaryVecPair buildBoundaries(const ConditionVecPair& leftRightBC, int spatialOrder, int nd);

	// Structure linking methods.
	static GridJoin* link(const int dim, GridCell* lcptr, GridCell* rcptr);
	static GridJoin* weakLink(const int dim, GridCell* lcptr, GridCell* rcptr);
	static void boundaryLink(Boundary* bptr, GridCell* const firstCell, const Coords& coreCells, std::array<JoinContainer, 3>& joins);

	// Geometry calculation methods.
	static double vol_cell(double rc, const Vec3& dx, Geometry geometry, int nd);
	static double area_join(const Vec3& xj, const int dim, const Vec3& dx, Geometry geometry, int nd);
};



#endif /* GRIDBUILDER_HPP_ */
