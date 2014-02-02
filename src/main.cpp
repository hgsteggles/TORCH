/* main.cpp */

#include <stdio.h>
#include <stdlib.h>
#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <cmath>
#include <limits>
#include "main.hpp"
void setupParameters(GridParameters& gpar, RadiationParameters& rpar, HydroParameters& hpar, PrintParameters& ppar, Scalings& spar);
void initialConditions(Grid3D* gptr, HydroDynamics& hydro, Radiation& rad, Scalings& scale, MPIHandler& mpih);
void march(double tmax, Grid3D* gptr, const Radiation& rad, const HydroDynamics& hydro, const InputOutput& io, MPIHandler& mpihandler);
double fluidStep(Grid3D* gptr, const Radiation& rad, const HydroDynamics& hydro, MPIHandler& mpih);
void deleteFileContents(string myString);
//void printParams(double runtime, Parameters& p);
bool SWEEPX = true;

int main (){
	/* DECLARE VARIABLES AND SET UP PARAMETERS */
	MPIHandler mpihandler;
	if (mpihandler.getRank() == 0)
		deleteFileContents("tmp/");
	GridParameters gpar;
	RadiationParameters rpar;
	HydroParameters hpar;
	PrintParameters ppar;
	Scalings scale;
	setupParameters(gpar, rpar, hpar, ppar, scale);
	Grid3D *gptr = new Grid3D(gpar, mpihandler);
	HydroDynamics hydro(hpar);
	Radiation rad(rpar);
	InputOutput io(ppar, scale);
	initialConditions(gptr, hydro, rad, scale, mpihandler);
	cout << "MAIN: initial conditions set. Starting march... " << '\n';
	double n_H = rpar.NHI / (1.0/(scale.L*scale.L*scale.L));
	//double trec = 1.0/(ALPHA*n_H);
	double RSinf = pow((3.0*rpar.SOURCE_S)/(4.0*PI*n_H*n_H*rpar.ALPHA_B), 1.0/3.0);
	double cII = sqrt(GAS_CONST*2.0*rpar.TMAX*(scale.P/scale.RHO)) / scale.V;
	double t_s = RSinf/cII;
	double tmax = 2*0.0001*14.0*t_s;
	tmax = 0.2;
	mpihandler.barrier();
	march(tmax, gptr, rad, hydro, io, mpihandler);
	delete gptr;
	/* PRINT PARAMS */
	//printParams(runTime, p);
	return 0;
}

void setupParameters(GridParameters& gpar, RadiationParameters& rpar, HydroParameters& hpar, PrintParameters& ppar, Scalings& scale){
	scale.set_LMT(14.6*PC2CM, 1.0*MO2G, 10000*YR2S); // set length, mass & time scalings.
	gpar.ND = 1; // No. of spatial dimensions.
	gpar.NU = 6; // No. of conserved variables involved in riemann solver.
	gpar.NR = 5; // No. of radiation variables not involved in riemann solver.
	gpar.NCELLS[0] = 1280; // No. of cells along 1st dimension.
	gpar.NCELLS[1] = 1; // No. of cells along 2nd dimension.
	gpar.NCELLS[2] = 1; // No. of cells along 3rd dimension.
	gpar.GEOMETRY = SPHERICAL; // Geometry of grid {CARTESIAN, CYLINDRICAL, SPHERICAL}.
	gpar.ORDER_S = 1; // Order of reconstruction in cell {0=CONSTANT, 1=LINEAR}.
	gpar.ORDER_T = 2; // Order of accuracy in a time step {1=FirstOrder, 2=SecondOrder}.
	rpar.K1 = 0.2; // [Mackey 2012 (table 1)] timestep constant.
	rpar.K2 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	rpar.K3 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	rpar.K4 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	rpar.P_I_CROSS_SECTION = 6.3E-18 / (scale.L*scale.L); // P.I. cross-sect (cm^2/scale).
	rpar.ALPHA_B = 2.59E-13 / (scale.L*scale.L*scale.L/scale.T); // Case B radiative recombination rate (cm^3 t^-1/scale).
	rpar.TAU_0 = 0.6; // Min. tau in nearest neighbour weights [Mellema et. al. 2006 (A.5)].
	rpar.SOURCE_S = 1e51 / (1.0/scale.T); // Source photon Luminosity (s^-1/scale).
	rpar.NHI = 4000; // Initial number density of neutral hydrogen (cm^-3).
	rpar.TMIN = 100 / (scale.P/scale.RHO); // Temperature fix for fully neutral gas.
	rpar.TMAX = 10000 / (scale.P/scale.RHO); // Temperature fix for fully ionized gas.
	rpar.SCHEME = IMPLICIT; // Ionization fraction integration scheme.
	rpar.H_MASS = H_MASS_G / scale.M; // Mass of hydrogen.
	hpar.GAMMA = 1.4; // Adiabatic gas constant.
	hpar.DFLOOR = 0.00000001; // Density floor to prevent negative density.
	hpar.PFLOOR = 0.00000001; // Pressure floor to prevent negative pressure.
	hpar.DTMAX = 0.001; // Maximum timestep.
	ppar.DIR_2D = "tmp/"; // Print directory for 2D grid data.
	ppar.DIR_IF = "tmp/"; // Print directory for ionization front data.
	ppar.PRINT2D_ON = true;
	ppar.PRINTIF_ON = false;
}

void initialConditions(Grid3D* gptr, HydroDynamics& hydro, Radiation& rad, Scalings& scale, MPIHandler& mpih){
	//Star star(0, 0, 0);
	//rad.addStar(star, gptr, mpih);
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
		/** SPITZER
		double rho = rad.NHI*H_MASS_G;
		double pre = (GAS_CONST/1.0)*rho*rad.TMIN*(scale.P/scale.RHO);
		cptr->Q[ihii] = 0.0;
		if ((cptr->xc[0]-1)*gptr->dx[0]*scale.L*CM2PC <  1.24989) {
			cptr->Q[ihii] = 1.0;
			pre = (GAS_CONST/0.5)*rho*rad.TMAX*(scale.P/scale.RHO);
		}
		cptr->Q[iden] = rho/scale.RHO;
		cptr->Q[ipre] = pre/scale.P;
		if(cptr == gptr->fsrc)
			cptr->Q[ihii] = 1.0;
		*/
		/** SOD */
		if (cptr->xc[0] < gptr->TOTNCELLS[0]/2) {
			cptr->Q[iden] = 1.0;
			cptr->Q[ipre] = 1.0;
		}
		else {
			cptr->Q[iden] = 0.125;
			cptr->Q[ipre] = 0.1;
		}
		hydro.UfromQ(cptr->U, cptr->Q);
	}
	if(rad.nstars > 0)
		rad.rayTrace(gptr);
	else
		cout << "NO SOURCES" << endl;
}
/**
 * Solves the state of radiation and hydrodynamic variables until the time tmax.
 * Marches the conserved variables, Grid3D::U, in a Grid3D object.
 * @param tmax The time (unitless) that march stops marching.
 * @param gptr A pointer to a Grid3D object, the object which is to be marched.
 * @param rad A reference to a Radiation object for solving the radiation field.
 * @param hydro A reference to a HydroDynamics object for solving the Riemann problem.
 * @param io A reference to an InputOutput object for printing.
 * @param mpih A reference to a HydroDynamics object.
 */
void march(double tmax, Grid3D* gptr, const Radiation& rad, const HydroDynamics& hydro, const InputOutput& io, MPIHandler& mpih){
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed;
	double dt = 0, t = 0;
	start = std::chrono::system_clock::now();
	cout << mpih.cname() << "MARCH: Marching solution...   0%" << endl;
	io.print2D(0, t, dt, gptr, mpih);
	if(rad.nstars > 0)
		io.printIF(t, rad, gptr, mpih);
	for(int step = 0, prc = 0; t < tmax; t += dt, step++){
		int prcnow = (int) (100*t/tmax);
		if (prcnow - prc >= 5){
			io.print2D(prcnow, t, dt, gptr, mpih);
			if(rad.nstars > 0)
				io.printIF(t, rad, gptr, mpih);
			cout << mpih.cname() << "MARCH: Marching solution...  ";
			if (prcnow == 5)
				cout << " ";
			cout << prcnow << "% " << dt << endl;
			prc = prcnow;
		}
		dt = fluidStep(gptr, rad, hydro, mpih);
		//cerr << "percent_done = " << t/tmax << endl;
	}
	io.print2D(100, t, dt, gptr, mpih);
	cout << mpih.cname() << "MARCH: Marching solution... " << 100 << "% " << dt << endl;
	end = std::chrono::system_clock::now();
	elapsed = end-start;
	mpih.barrier();
	cout << "MARCH: Took " << elapsed.count() << " seconds." << endl;
}
double fluidStep(Grid3D* gptr, const Radiation& rad, const HydroDynamics& hydro, MPIHandler& mpih){
	double dt_hydro, dt_rad, dt, dth, IF = 0;
	hydro.globalQfromU(gptr);
	hydro.updateBoundaries(gptr);
	if(gptr->ORDER_S == 1)
		hydro.reconstruct(gptr);
	hydro.calcFluxes(gptr);
	dt_hydro = hydro.CFL(gptr);
	dt_rad = rad.getTimeStep(dt_hydro, gptr);
	dt = min(dt_hydro, dt_rad);
	if (mpih.nProcessors() > 1)
		mpih.minimum(dt);
	if(gptr->ORDER_T == 2){
		dth = dt/2.0;
		hydro.globalWfromU(gptr);
	}
	else
		dth = dt;
	if (mpih.nProcessors() == 1)
		rad.transferRadiation(dth, IF, gptr);
	else
		rad.transferRadiation(dth, IF, gptr, mpih);
	hydro.applySrcTerms(dth, gptr, rad);
	hydro.advSolution(dth, gptr);
	hydro.fixSolution(gptr);
	if(gptr->ORDER_T == 2){
		hydro.globalQfromU(gptr);
		hydro.updateBoundaries(gptr);
		if(gptr->ORDER_S == 1)
			hydro.reconstruct(gptr);
		hydro.calcFluxes(gptr);
		hydro.globalUfromW(gptr);
		if (mpih.nProcessors() == 1)
			rad.transferRadiation(dt, IF, gptr);
		else
			rad.transferRadiation(dt, IF, gptr, mpih);
		hydro.applySrcTerms(dt, gptr, rad);
		hydro.advSolution(dt, gptr);
		hydro.fixSolution(gptr);
	}
	return dt;
}
////// deleteFileContents deletes all files within a directory
void deleteFileContents(string folder){
	struct dirent *next_file;
	DIR *dir; // These are data types defined in the "dirent" header.
	char filepath[256];
	dir = opendir(folder.c_str() );
	if(dir != NULL){
		while((next_file = readdir(dir) )){
			// build the full path for each file in the folder
			sprintf(filepath, "%s/%s", folder.c_str(), next_file->d_name);
			remove(filepath);
		}
	}
	cout << "IO: deleted all files in " << folder << '\n';
}
