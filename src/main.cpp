/* main.cpp */

#include "parameters.hpp"
#include "integrator.hpp"
#include "grid3d.hpp"
#include "gridcell.hpp"
#include "hydro.hpp"
#include "io.hpp"
#include "rtmodule.hpp"
#include "mpihandler.hpp"
#include "external.hpp"

//#include <google/profiler.h> //ProfilerStart, ProfilerStop
#include <dirent.h>

void setupParameters(GridParameters& gpar, RadiationParameters& rpar, HydroParameters& hpar, PrintParameters& ppar, Scalings& spar);
void deleteFileContents(const std::string& myString);
bool SWEEPX = true;

int main (){
	/* SIZES
	cout << "Size of GridCell = " << sizeof(GridCell) << '\n';
	cout << "Size of Grid3D = " << sizeof(Grid3D) << '\n';
	cout << "Size of Boundary = " << sizeof(Boundary) << '\n';
	cout << "Size of ExternalBoundary = " << sizeof(ExternalBoundary) << '\n';
	cout << "Size of HydroDynamics = " << sizeof(HydroDynamics) << '\n';
	cout << "Size of InputOutput = " << sizeof(InputOutput) << '\n';
	cout << "Size of MPIHandler = " << sizeof(MPIHandler) << '\n';
	cout << "Size of Partition = " << sizeof(Partition) << '\n';
	cout << "Size of Star = " << sizeof(Star) << '\n';
	exit(EXIT_FAILURE);
	*/
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
	Integrator integrator(gpar, rpar, hpar, ppar, scale, mpihandler);
	if (mpihandler.getRank() == 0)
		std::cout << "INTEGRATOR: initialising..." << '\n';
	integrator.init();
	integrator.io->addPrintTime(1);
	integrator.io->addPrintTime(2);
	integrator.io->addPrintTime(4);
	integrator.io->addPrintTime(8);
	integrator.io->addPrintTime(16);
	integrator.io->addPrintTime(28);
	InputOutput::debugging = false;

	double n_H = rpar.NHI / (1.0/(scale.L*scale.L*scale.L));
	//double trec = 1.0/(ALPHA*n_H);
	double RSinf = pow((3.0*rpar.SOURCE_S)/(4.0*PI*n_H*n_H*rpar.ALPHA_B), 1.0/3.0);
	double cII = sqrt(GAS_CONST*2.0*rpar.TMAX*(scale.P/scale.RHO)) / scale.V;
	double t_s = RSinf/cII;
	double tmax = 2*0.0001*14.0*t_s;
	tmax = 28.2;
	mpihandler.barrier();
	if (mpihandler.getRank() == 0)
		std::cout << "INTEGRATOR: starting march..." << '\n';
	//ProfilerStart("./main.prof");
	integrator.march(tmax);
	//ProfilerStop();
	/* PRINT PARAMS */
	//printParams(runTime, p);
	return 0;
}

void setupParameters(GridParameters& gpar, RadiationParameters& rpar, HydroParameters& hpar, PrintParameters& ppar, Scalings& scale){
	scale.set_LMT(1.256*PC2CM, 1.0*MO2G, 5000*YR2S); // set length, mass & time scalings.
	gpar.ND = 3; // No. of spatial dimensions.
	gpar.NU = 6; // No. of conserved variables involved in riemann solver.
	gpar.NR = 5; // No. of radiation variables not involved in riemann solver.
	gpar.NCELLS[0] = 128; // No. of cells along 1st dimension.
	gpar.NCELLS[1] = 128; // No. of cells along 2nd dimension.
	gpar.NCELLS[2] = 128; // No. of cells along 3rd dimension.
	gpar.GEOMETRY = CARTESIAN; // Geometry of grid {CARTESIAN, CYLINDRICAL, SPHERICAL}.
	gpar.ORDER_S = 1; // Order of reconstruction in cell {0=CONSTANT, 1=LINEAR}.
	gpar.ORDER_T = 2; // Order of accuracy in a time step {1=FirstOrder, 2=SecondOrder}.
	gpar.LBCondition[0] = REFLECTING;
	gpar.LBCondition[1] = REFLECTING;
	gpar.LBCondition[2] = REFLECTING;
	gpar.RBCondition[0] = REFLECTING;
	gpar.RBCondition[1] = REFLECTING;
	gpar.RBCondition[2] = REFLECTING;
	rpar.K1 = 0.02; // [Mackey 2012 (table 1)] timestep constant.
	rpar.K2 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	rpar.K3 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	rpar.K4 = 0.0; // [Mackey 2012 (table 1)] timestep constant.
	rpar.P_I_CROSS_SECTION = 6.3E-18 / (scale.L*scale.L); // P.I. cross-sect (cm^2/scale).
	rpar.ALPHA_B = 2.7e-13 / (scale.L*scale.L*scale.L/scale.T); // Case B radiative recombination rate (cm^3 t^-1/scale).
	rpar.TAU_0 = 0.6; // Min. tau in nearest neighbour weights [Mellema et. al. 2006 (A.5)].
	rpar.SOURCE_S = 1e49 / (1.0/scale.T); // Source photon Luminosity (s^-1/scale).
	rpar.NHI = 5.21e-21/H_MASS_G; // Initial number density of neutral hydrogen (cm^-3).
	rpar.TMIN = 100 / (scale.P/scale.RHO); // Temperature fix for fully neutral gas.
	rpar.TMAX = 10000 / (scale.P/scale.RHO); // Temperature fix for fully ionized gas.
	rpar.SCHEME = IMPLICIT; // Ionization fraction integration scheme.
	rpar.H_MASS = H_MASS_G / scale.M; // Mass of hydrogen.
	hpar.GAMMA = 1.0001; // Adiabatic gas constant.
	hpar.DFLOOR = 0.00000001; // Density floor to prevent negative density.
	hpar.PFLOOR = 0.00000001; // Pressure floor to prevent negative pressure.
	hpar.DTMAX = 100; // Maximum timestep.
	ppar.DIR_2D = "tmp/"; // Print directory for 2D grid data.
	ppar.DIR_IF = "tmp/"; // Print directory for ionization front data.
	ppar.PRINT2D_ON = true;
	ppar.PRINTIF_ON = false;
}

////// deleteFileContents deletes all files within a directory
void deleteFileContents(const std::string& folder){
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
	std::cout << "IO: deleted all files in " << folder << '\n';
}
