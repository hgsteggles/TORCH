/**
 * Provides the Radiation class.
 * @file radiation.hpp
 * @author Harrison Steggles
 * @date 29/01/2014 - The first version.
 * @date 30/01/2014 - Star passed into addStar can now have any location. Pointer to first GridCell in causal
 * loop (fcausal) is nearest to Star with preference to positive directions.
 * @date 31/01/2014 - Star can only snap to grid vertex at the corner of a GridCell nearest the origin.
 * Coordinates passed in are the grid coordinates of the GridCell. Therefore coordinate transforms should
 * be considered when comparing locations. GridCell::xc on Star x axis is = GridCell::xc + 0.5.
 * @date 31/01/2014 - Fixed bug in first if statement in Radiation::updateTauSC. Now escapes method if Star
 * is located in one of GridCells corners.
 * @date 04/01/2014, arguments now passed by const reference when appropriate.
 */

#ifndef RTMODULE_H
#define RTMODULE_H

#include <stdlib.h>
#include <vector>

class GridCell;
class Grid3D;
class RadiationParameters;
class Star;
class MPIHandler;
class Partition;

/**
 * @class Radiation
 * @brief Contains parameters and methods for calculating the radiation field around a Star object in a grid.
 * @version 0.4, 04/02/2014
 */
class Radiation{
public:
	Grid3D* gptr;
	double K1, K2, K3, K4, P_I_CROSS_SECTION, ALPHA_B, TAU_0, SOURCE_S, NHI, TMIN, TMAX, SCHEME, H_MASS;
	std::vector<Star> stars;
	int nstars;
	Radiation(const RadiationParameters& rp, Grid3D* grid);
	void addStar(Star src, MPIHandler& mpih);
	int getRayPlane(GridCell* cptr, const int& starID) const;
	void updateTauSC(bool average, GridCell* cptr, const int& starID) const;
	double temperature(GridCell* cptr) const;
	double alphaB(GridCell* cptr) const;
	double cellPathLength(GridCell* cptr, const int& starID) const;
	double shellVolume(GridCell* cptr, const int& starID) const;
	double PIrate(const double& frac, const double& T, const double& dT, GridCell* cptr, const int& starID) const;
	double HIIfracDot(const double& A_pi, const double& frac, GridCell* cptr) const;
	void update_dtau(GridCell* cptr, const int& starID) const;
	double radHeatCool(const double& dt, GridCell* cptr) const;
	void doric(const double& dt, double& frac, double& frac_av, const double& A_pi, GridCell* cptr) const;
	void update_HIIfrac(const double& dt, GridCell* cptr, const int& starID) const;
	double getTimeStep(const double& dt_dyn) const;
	void transferRadiation(const double& dt, double& IF) const;
	void transferRadiation(const double& dt, double& IF, MPIHandler& mpih) const;
	void sendColumnDensities(MPIHandler& mpih) const;
	void addSource(const double& x, const double& y, const double& z) const;
	void printIF(const int& starID, const double& t) const;
	void rayTrace() const;
	GridCell* locateStar(Star star) const;
	bool isStar(GridCell* cptr) const;
};

#endif
