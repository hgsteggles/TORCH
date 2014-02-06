/**
 * Provides the GridCell and GridJoin classes.
 * @file gridcell.hpp
 *
 * @author Harrison Steggles
 * @date 13/01/2014, the first version.
 * @date 16/01/2014, removed old Boundary class.
 * @date 04/01/2014, GridCell now has pointer to next GridCell in causal iteration (nextcausal). UL & UR unused states removed.
 * @date 04/01/2014, arguments now passed by const reference when appropriate.
 */

#ifndef GRIDCELL_H
#define GRIDCELL_H

#include "constants.hpp"

class GridJoin;
class GhostCell;
/**
 * @class GridCell
 * @brief The GridCell class holds the fluid state information in a single grid cell of Grid3D.
 * A GridCell object holds links (for fast traversal) to nearest neighbours and GridJoin objects that lie between it and its neighbours.
 * @version 0.3, 28/01/2014
 */
class GridCell{
public:
	/**
	 * @brief Default GridCell constructor.
	 * Provides all attributes with safe values.
	 */
	GridCell();
	/**
	 * @brief GridCell constructor.
	 * Does the same as the default constructor except that the grid coordinates are passed.
	 * @param i The ith cooordinate of the GridCell in a Grid3D.
	 * @param j The jth cooordinate of the GridCell in a Grid3D.
	 * @param k The kth cooordinate of the GridCell in a Grid3D.
	 */
	GridCell(const int& i, const int& j, const int& k);
	/**
	 * @brief GridCell destructor.
	 * Does NOT delete objects that this GridCell object points to.
	 */
	~GridCell();
	void printInfo();
	/**
	 * @brief Setter for GridCell::U.
	 * @param index
	 * @param value
	 */
	void set_U(const int& index, const double& value);
	/**
	 * @brief Setter for GridCell:xc.
	 * @param x The x grid coordinate.
	 * @param y The y grid coordinate.
	 * @param z The z grid coordinate.
	 */
	void set_xcs(const int& x, const int& y, const int& z);
	/**
	 * @brief Getter for GridCell::xc.
	 * @param i The grid coordinate to be returned.
	 * @return The location of this GridCell object on grid coordinate i.
	 */
	int get_xc(const int& i);
	/**
	 * @brief Getter for GridCell::U.
	 * @param index The index for the fluid variable to be returned.
	 * @return The value of the fluid variable.
	 */
	double get_U(const int& index);
	/**
	 * @brief Returns the temperature of this GridCell object.
	 * @return Temperature.
	 */
	double temperature();

	static int s_total; //!< A count of all GridCell objects created in the program.
	GridJoin* rjoin[3]; //!< Contains pointers to GridJoin objects that lie on the right side of this GridCell.
	GridJoin* ljoin[3]; //!< Contains pointers to GridJoin objects that lie on the left side of this GridCell.
	GridCell* right[3]; //!< Contains pointers to GridCell objects that lie on the right side of this GridCell.
	GridCell* left[3]; //!< Contains pointers to GridCell objects that lie on the left side of this GridCell.
	GridCell* next; //!< Points to the next GridCell object in a list that Grid3D uses to clean up GridCell objects.
	GridCell* nextcausal; //!< Points to the next GridCell object in a causal list that starts at a source of radiation.
	double U[NU]; //!< Contains conservative fluid variable values.
	double Q[NU]; //!< Contains primitive fluid variable values.
	double W[NU]; //!< Contains a copy of GridCell::U for 2nd order time-stepping.
	double R[NR]; //!< Contains radiation variable values: optical depth in cell and along path of the ray from source.
	int xc[3]; //!< Grid coordinates for this GridCell.
	double QL[3][NU]; //!< Reconstructed states on left faces.
	double QR[3][NU]; //!< Reconstructed states on right faces.
	double vol; //!< Volume of GridCell.
};

/**
 * @class GridJoin
 * @brief The GridJoin class contains flux information.
 *
 * A GridJoin object holds links to GridCell objects that lie either side of it and contains flux information. The GridCell objects either side
 * of this will receive the flux in this GridJoin which is calculated by the HLLC Riemann Solver in a HydroDynamics object.
 * @version 0.3, 28/01/2014
 */
class GridJoin{
public:
	/**
	 * @brief The default GridJoin constructor.
	 * Provides all attributes with safe values.
	 */
	GridJoin();
	/**
	 * @brief The Gridjoin destructor.
	 */
	~GridJoin();
	static int s_total; //!< A count of all GridJoin objects created in the program.
	GridCell* lcell; //!< Pointer to GridCell object on the left.
	GridCell* rcell; //!< Pointer to GridCell object on the right.
	GridJoin* next; //!< Points to the next GridJoin object in a list that Grid3D uses to clean up GridJoin objects.
	double F[NU]; //!< Contains flux to be added to GridJoin::rcell and subtracted from GridJoin::lcell fluid variables (GridCell::U).
	double xj[NU]; //!< The grid coordinates of the GridJoin object in a Grid3D object.
	double area; //!< The area of the GridJoin.
};

#endif
