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
 * @date 04/02/2014 - arguments now passed by const reference when appropriate.
 * @date 13/02/2014 - 28% speedup achieved by adding GridCell member variables that hold cell path length,
 * pointers to nearest neighbours and the weights associated with them instead of calculating them at each
 * time step.
 * @date 11/03/2014 - bug in transferRadiation(const double&, double&, MPIHandler&). Was accessing boundaries[3]
 * instead of boundaries[1].
 * @date 14/04/2014 - Radiation parameters now all held within a RadiationParameters object.
 * @date 14/04/2014 - function updateColDen added that calculates column density for total hydrogen density
 * (rather than ionized hydrogen as in updateTauSC).
 * @date 14/04/2014 - function collIonRate added in anticipation of full support for collisional ionizations.
 * @date 01/05/2014 - getTimeStep is now calcTimeStep. More accurate named after what it does.
 * @date 28/05/2014 - new methods added to update source terms instead of directly modifying the conservative fluid variables.
 */

#ifndef RADIATION_H
#define RADIATION_H

#include <stdlib.h>
#include <vector>

class GridCell;
class Grid3D;
class RadiationParameters;
class HydroParameters;
class Scalings;
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
	RadiationParameters* rparams;
	Scalings* scale;
	Grid3D* gptr;
	std::vector<Star> stars;
	int nstars;

	Radiation(const RadiationParameters& rp, const Scalings& sc, Grid3D* grid);

	//Initialization methods.
	int getRayPlane(GridCell* cptr, const int& starID) const;
	void calculateNearestNeighbours(const int& starID) const;
	void rayTrace() const;
	double cellPathLength(GridCell* cptr, const int& starID) const;
	double shellVolume(GridCell* cptr, const int& starID) const;


	//Calculation methods.
	void doric(const double& dt, double& frac, double& frac_av, const double& A_pi, GridCell* cptr) const;
	double temperature(GridCell* cptr) const;
	double alphaB(const double& T) const;
	double collIonRate(const double& Temperature);
	double PIrate(const double& frac, const double& T, const double& dT, GridCell* cptr, const int& starID) const;
	double HIIfracDot(const double& A_pi, const double& frac, GridCell* cptr) const;
	double calcTimeStep(const double& dt_dyn) const;

	//Update methods.
	void updateTauSC(bool average, GridCell* cptr, const int& starID) const;
	void updateTauSC2(bool average, GridCell* cptr, const double& dist2) const;
	void updateColDen(GridCell* cptr, const double& dist2) const;
	void update_dtau(GridCell* cptr, const int& starID) const;
	void update_HIIfrac(const double& dt, GridCell* cptr, const int& starID) const;

	//Integration methods.
	void transferRadiation(const double& dt) const;
	void transferRadiation(const double& dt, MPIHandler& mpih) const;
	void applySrcTerms(const double& dt, const HydroParameters* hparams);
	void updateSrcTerms(const double& dt, const HydroParameters* hparams);

	//Misc. methods.
	void addStar(Star src, MPIHandler& mpih, bool snapToFace[6]);
	bool isStar(GridCell* cptr) const;
};

#endif
