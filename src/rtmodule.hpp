/* rtmodule.h */

#ifndef RTMODULE_H
#define RTMODULE_H

#include "gridcell.hpp"
#include "grid3d.hpp"
#include "parameters.hpp"
#include <stdlib.h>

class Radiation{
public:
	double K1, K2, K3, K4, P_I_CROSS_SECTION, ALPHA_B, TAU_0, SOURCE_S, NHI, TMIN, TMAX, SCHEME, H_MASS;
	Radiation(const RadiationParameters& rp);
	int getRayPlane(GridCell* cptr, GridCell* srcptr, Grid3D* gptr) const;
	void updateTauSC(bool average, GridCell* cptr, GridCell* srcptr, Grid3D* gptr) const;
	double temperature(GridCell* cptr) const;
	double alphaB(GridCell* cptr) const;
	double cellPathLength(GridCell* cptr, GridCell* srcptr, Grid3D* gptr) const;
	double shellVolume(GridCell* cptr, GridCell* srcptr, Grid3D* gptr) const;
	double PIrate(double frac, double T, double dT, GridCell* cptr, GridCell* srcptr, Grid3D* gptr) const;
	double HIIfracDot(double A_pi, double frac, GridCell* cptr) const;
	void update_dtau(GridCell* cptr, GridCell* srcptr, Grid3D* gptr) const;
	double radHeatCool(double dt, GridCell* cptr) const;
	void doric(double dt, double& frac, double& frac_av, double A_pi, GridCell* cptr) const;
	void update_HIIfrac(double dt, GridCell* cptr, GridCell* srcptr, Grid3D* gptr) const;
	void update_HIIfrac(double dt, GridCell* cptr, GridCell* srcptr) const;
	double getTimeStep(double dt_dyn, GridCell* srcptr, Grid3D* gptr) const;
	void transferRadiation(double dt, double& IF, Grid3D* gptr) const;
	void addSource(double x, double y, double z) const;
	void printIF(GridCell* srcptr, Grid3D* gptr, double t) const;
	void rayTrace(Grid3D* gptr) const;
};

#endif
