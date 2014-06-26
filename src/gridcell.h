/**
 * Provides the GridCell and GridJoin classes.
 * @file gridcell.h
 *
 * @author Harrison Steggles
 * @date 13/01/2014, the first version.
 * @date 16/01/2014, removed old Boundary class.
 * @date 04/02/2014, GridCell now has pointer to next GridCell in causal iteration (nextcausal). UL & UR unused states removed.
 * @date 04/02/2014, arguments now passed by const reference when appropriate.
 * @date 12/02/2014, printInfo method added for debugging.
 * @date 13/02/2014 - GridCell now has pointers to nearest neighbouring GridCells, their associated weights for interpolating
 * column density, the shell volume and cell crossing length used in radiative transfer methods. Get a speed-up of 27% by storing
 * these rather than calculating them at each step.
 */

#ifndef GRIDCELL_H
#define GRIDCELL_H

#include "constants.h"

class GridJoin;
class GhostCell;
/**
 * @class GridCell
 *
 * @brief The GridCell class holds the fluid state information in a single grid cell of Grid3D.
 *
 * A GridCell object holds links (for fast traversal) to nearest neighbours and GridJoin objects that lie between it and its neighbours.
 *
 * @version 0.7, 13/06/2014
 */
class GridCell{
public:
	static int s_total; //!< A count of all GridCell objects created in the program.
	GridJoin* rjoin[3]; //!< Contains pointers to GridJoin objects that lie on the right side of this GridCell.
	GridJoin* ljoin[3]; //!< Contains pointers to GridJoin objects that lie on the left side of this GridCell.
	GridCell* right[3]; //!< Contains pointers to GridCell objects that lie on the right side of this GridCell.
	GridCell* left[3]; //!< Contains pointers to GridCell objects that lie on the left side of this GridCell.
	GridCell* next; //!< Points to the next GridCell object in a list that Grid3D uses to clean up GridCell objects.
	GridCell* nextcausal; //!< Points to the next GridCell object in a causal list that starts at a source of radiation.
	double UDOT[NU]; //!< Contains rate of change of conservative fluid variable values.
	double U[NU]; //!< Contains conservative fluid variable values.
	double Q[NU]; //!< Contains primitive fluid variable values.
	double W[NU]; //!< Contains a copy of GridCell::U for 2nd order time-stepping.
	double R[NR]; //!< Contains radiation variable values: optical depth in cell and along path of the ray from source.
	int xc[3]; //!< Grid coordinates for this GridCell.
	double QL[3][NU]; //!< Reconstructed states on left faces.
	double QR[3][NU]; //!< Reconstructed states on right faces.
	double vol; //!< Volume of GridCell.
	double ds;
	double shellVol;
	GridCell* NN[4];
	double NN_weights[4];

	//Structors.
	GridCell();
	~GridCell();

	//Debugging.
	void printInfo();

	//Setters and Getters.
	void set_U(const int& index, const double& value);
	void set_xcs(const int& x, const int& y, const int& z);
	int get_xc(const int& i);
	double get_U(const int& index);

	//Misc. methods.
	double temperature();
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
	static int s_total; //!< A count of all GridJoin objects created in the program.
	GridCell* lcell; //!< Pointer to GridCell object on the left.
	GridCell* rcell; //!< Pointer to GridCell object on the right.
	GridJoin* next; //!< Points to the next GridJoin object in a list that Grid3D uses to clean up GridJoin objects.
	double F[NU]; //!< Contains flux to be added to GridJoin::rcell and subtracted from GridJoin::lcell fluid variables (GridCell::U).
	double xj[NU]; //!< The grid coordinates of the GridJoin object in a Grid3D object.
	double area; //!< The area of the GridJoin.

	GridJoin();
	~GridJoin();
};

#endif
