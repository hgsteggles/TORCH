/**
 * Provides the InputOutput class.
 * @file io.hpp
 *
 * @author Harrison Steggles
 * @date 13/01/2014, the first version.
 */

#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <map>
#include "parameters.hpp"
#include "grid3d.hpp"
#include "rtmodule.hpp"
#include "mpihandler.hpp"

/**
 * @class InputOutput
 * @brief Handler for all input/output operations.
 * @version 0.3, 28/01/2014
 */
class InputOutput{
public:
	string DIR_2D; //!< Directory the output of InputOutput::print2D is aimed.
	string DIR_IF; //!< Directory the output of InputOutput::printIF is aimed.
	bool PRINT2D_ON; //!< Toggles whether InputOutput::print2D will print anything.
	bool PRINTIF_ON; //!< Toggles whether InputOutput::printIF will print anything.
	double L_SCALE; //!< Length scale.
	double M_SCALE; //!< Mass scale.
	double T_SCALE; //!< Time scale.
	double V_SCALE; //!< Velocity scale.
	double RHO_SCALE; //!< Density scale.
	double P_SCALE; //!< Pressure scale.
	double E_SCALE; //!< Energy scale.
	/**
	 * @brief InputOutput constructor.
	 * @param rp Parameters to pass in.
	 * @param sc Scalings to pass in for printing physical units.
	 */
	InputOutput(const PrintParameters& rp, const Scalings& sc);
	/**
	 * @brief Prints primitive variables for a 2D slice of the simulation grid in the Grid3D object pointed to by gptr.s
	 * @param step Number to append to the output filename.
	 * @param t Current simulation time.
	 * @param dt Current delta-time.
	 * @param gptr Pointer to Grid3D object to print.
	 * @param mpih Provides MPI information for printing with multiple cores.
	 */
	void print2D(int step, double t, double dt, Grid3D* gptr, MPIHandler& mpih) const;
	/**
	 * @brief Calculates and prints the location of the ionization front from a Star that is located at the origin.
	 * @param t Current simulation time.
	 * @param rad Radiation provides this method with radiation parameters.
	 * @param gptr The grid on which the Star is located.
	 * @param mpih Provides MPI information for calculating and printing with multiple cores.
	 */
	void printIF(double t, const Radiation& rad, Grid3D* gptr, MPIHandler& mpih) const;
	/**
	 * @brief Prints all parameters associated with this simulation. UNFINISHED.
	 * @param runtime Simulation time.
	 * @param gpar Grid3D parameters for printing.
	 * @param rpar Radiation parameters for printing.
	 * @param hpar Hydrodynamics parameters for printing.
	 * @param ppar I/O parameters for printing.
	 * @param scale Scales for printing.
	 */
	void printParams(double runtime, const GridParameters& gpar, const RadiationParameters& rpar, const HydroParameters& hpar, const PrintParameters& ppar, const Scalings& scale) const;
	void fileToMap(const string& myString, std::map<double,double> myMap) const;
};


#endif
