/**
 * @file main.cpp
 */

#include "parameters.h"
#include "parameterparser.h"
#include "integrator.h"
#include "grid3d.h"
#include "gridcell.h"
#include "hydro.h"
#include "io.h"
#include "radiation.h"
#include "mpihandler.h"
#include "external.h"

//#include <google/profiler.h> //ProfilerStart, ProfilerStop
#include <dirent.h>
#include <cmath>

#define _TO_STRING_(A) #A

void setupParameters(IntegrationParameters& ipar, GridParameters& gpar, RadiationParameters& rpar, HydroParameters& hpar, PrintParameters& ppar, Scalings& spar);
void setupParametersOld(IntegrationParameters& ipar, GridParameters& gpar, RadiationParameters& rpar, HydroParameters& hpar, PrintParameters& ppar, Scalings& spar);
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
	IntegrationParameters ipar;
	GridParameters gpar;
	RadiationParameters rpar;
	HydroParameters hpar;
	PrintParameters ppar;
	Scalings scale;

	for (int i = 0; i < mpihandler.nProcessors(); ++i) {
		if (mpihandler.getRank() == i)
			setupParameters(ipar, gpar, rpar, hpar, ppar, scale);
		mpihandler.barrier();
	}

	Integrator integrator(ipar, gpar, rpar, hpar, ppar, scale, mpihandler);
	if (mpihandler.getRank() == 0)
		std::cout << "INTEGRATOR: initialising..." << '\n';
	integrator.init(scale);
	//integrator.io->addPrintTime(1);
	//integrator.io->addPrintTime(2);
	//integrator.io->addPrintTime(4);
	//integrator.io->addPrintTime(8);
	//integrator.io->addPrintTime(16);
	//integrator.io->addPrintTime(28);
	InputOutput::debugging = false;

	double photon_rate = 0;
	if (rpar.vStarParams.size() != 0)
		photon_rate = rpar.vStarParams[0].PHOTON_RATE;
	double n_H = scale.toCodeUnits(rpar.NHI, 0, -3, 0);
	//double trec = 1.0/(rpar.ALPHA_B*n_H);
	double RSinf = std::pow((3.0*photon_rate)/(4.0*PI*n_H*n_H*rpar.ALPHA_B), 1.0/3.0);
	double cII = rpar.SPEC_GAS_CONST*2.0*rpar.THII;
	double t_s = RSinf/cII;
	//double tmax = 4*16.0*t_s;
	double tmax = 20;//*30;
	//tmax = 2.0;
	//int nsteps = 40000;
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

void setupParameters(IntegrationParameters& ipar, GridParameters& gpar, RadiationParameters& rpar,
		HydroParameters& hpar, PrintParameters& ppar, Scalings& scale) {
	ParameterParser parser;
	parser.parseParameters("refdata/parameters.xml", ipar, gpar, hpar, rpar, ppar);
	scale.set_rho_pressure_time(5.21e-21, 1.0e-10, 5000*YR2S);

	//ipar.print();
	//gpar.print();
	//rpar.print();
	//hpar.print();
	//ppar.print();
	//scale.print();

	gpar.SIDE_LENGTH = scale.toCodeUnits(gpar.SIDE_LENGTH, 0, 1, 0);
	rpar.P_I_CROSS_SECTION = scale.toCodeUnits(rpar.P_I_CROSS_SECTION, 0, 2, 0); // P.I. cross-sect (cm^2/scale).
	rpar.ALPHA_B = scale.toCodeUnits(rpar.ALPHA_B, 0, 3, -1); // Case B radiative recombination rate (cm^3 t^-1/scale).
	rpar.NHI /= H_MASS_G;
	rpar.H_MASS = scale.toCodeUnits(H_MASS_G, 1, 0, 0); // Mass of hydrogen.
	rpar.SPEC_GAS_CONST = scale.toCodeUnits(8.314462e7, 0, 2, -2); // Specific gas constant [g-1.ergs.mol-1.K-1].
	for (int i = 0; i < rpar.vStarParams.size(); ++i) {
		StarParameters& spar = rpar.vStarParams[i];
		spar.PHOTON_ENERGY = scale.toCodeUnits(spar.PHOTON_ENERGY, 1, 2, -2); // Source photon energy (erg/scale).
		spar.PHOTON_RATE = scale.toCodeUnits(spar.PHOTON_RATE, 0, 0, -1); // Source photon Luminosity (s^-1/scale).
	}
}

void setupParametersOld(IntegrationParameters& ipar, GridParameters& gpar, RadiationParameters& rpar, HydroParameters& hpar, PrintParameters& ppar, Scalings& scale){
	ipar.ORDER_S = 1; // Order of reconstruction in cell {0=CONSTANT, 1=LINEAR}.
	ipar.ORDER_T = 2; // Order of accuracy in a time step {1=FirstOrder, 2=SecondOrder}.
	ipar.DT_MAX = 100; // Maximum timestep.

	//scale.set_LMT(1.0*MO2G, 1.256*PC2CM, 5000*YR2S); // set length, mass & time scalings.
	//scale.set_mass_length_time(10*MO2G, 1.256*PC2CM, 5000*YR2S); // set length, mass & time scalings.
	scale.set_rho_pressure_time(5.21e-21, 1.0e-10, 5000*YR2S);

	gpar.ND = 2; // No. of spatial dimensions.
	gpar.NCELLS[0] = 128; // No. of cells along 1st dimension.
	gpar.NCELLS[1] = 128; // No. of cells along 2nd dimension.
	gpar.NCELLS[2] = 1; // No. of cells along 3rd dimension.
	gpar.GEOMETRY = CYLINDRICAL; // Geometry of grid {CARTESIAN, CYLINDRICAL, SPHERICAL}.
	gpar.SIDE_LENGTH = scale.toCodeUnits(1.256*PC2CM, 0, 1, 0);
	gpar.LBCondition[0] = REFLECTING;
	gpar.LBCondition[1] = REFLECTING;
	gpar.LBCondition[2] = REFLECTING;
	gpar.RBCondition[0] = REFLECTING;
	gpar.RBCondition[1] = REFLECTING;
	gpar.RBCondition[2] = REFLECTING;

	rpar.K1 = 0.2; // [Mackey 2012 (table 1)] timestep constant for radiative transfer.
	rpar.K2 = 0.0; // [Mackey 2012 (table 1)] timestep constant for radiative transfer.
	rpar.K3 = 0.0; // [Mackey 2012 (table 1)] timestep constant for radiative transfer.
	rpar.K4 = 0.0; // [Mackey 2012 (table 1)] timestep constant for radiative transfer.
	rpar.P_I_CROSS_SECTION = scale.toCodeUnits(6.3E-18, 0, 2.0, 0); // P.I. cross-sect (cm^2/scale).
	rpar.ALPHA_B = scale.toCodeUnits(2.7e-13, 0, 3.0, -1); // Case B radiative recombination rate (cm^3 t^-1/scale).
	rpar.TAU_0 = 0.6; // Min. tau in nearest neighbour weights [Mellema et. al. 2006 (A.5)].
	rpar.NHI = 5.21e-21/H_MASS_G; // Initial number density of neutral hydrogen (cm^-3).
	rpar.THI = 100; // Temperature fix for fully neutral gas.
	rpar.THII = 10000; // Temperature fix for fully ionized gas.
	rpar.SCHEME = IMPLICIT; // Ionization fraction integration scheme.
	rpar.H_MASS = scale.toCodeUnits(H_MASS_G, 1, 0, 0); // Mass of hydrogen.
	rpar.SPEC_GAS_CONST = scale.toCodeUnits(8.314462e7, 0, 2, -2); // Specific gas constant [g-1.ergs.mol-1.K-1].

	hpar.GAMMA = 1.4; // Adiabatic gas constant.
	hpar.DFLOOR = 0.00000001; // Density floor to prevent negative density.
	hpar.PFLOOR = 0.00000001; // Pressure floor to prevent negative pressure.

	ppar.DIR_2D = "tmp/"; // Print directory for 2D grid data.
	ppar.DIR_IF = "tmp/"; // Print directory for ionization front data.
	ppar.PRINT2D_ON = true;
	ppar.PRINTIF_ON = true;
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
