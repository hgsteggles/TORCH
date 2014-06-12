/**
 * Provides the InputOutput class.
 * @file io.hpp
 *
 * @author Harrison Steggles
 * @date 13/01/2014 - the first version.
 * @date 04/01/2014 - arguments now passed by const reference when appropriate.
 * @date 12/02/2014 - progress bar added which outputs progress to command line and sets the frequency for
 * printing data, initProgressBar needs to be called first and progressBar() needs to be called at the end.
 */

#ifndef IO_H
#define IO_H

#include "parameters.hpp"
#include "grid3d.hpp"
#include "radiation.hpp"
#include "mpihandler.hpp"

#include <map>

class PrintParameters;
class Scalings;
class Grid3D;
class MPIHandler;
class Radiation;

/**
 * @class InputOutput
 * @brief Handler for all input/output operations.
 * @version 0.5, 24/02/2014
 */
class InputOutput{
public:
	std::vector<double> printTimes;
	std::string DIR_2D; //!< Directory the output of InputOutput::print2D is aimed.
	std::string DIR_IF; //!< Directory the output of InputOutput::printIF is aimed.
	bool PRINT2D_ON; //!< Toggles whether InputOutput::print2D will print anything.
	bool PRINTIF_ON; //!< Toggles whether InputOutput::printIF will print anything.
	double L_SCALE; //!< Length scale.
	double M_SCALE; //!< Mass scale.
	double T_SCALE; //!< Time scale.
	double V_SCALE; //!< Velocity scale.
	double RHO_SCALE; //!< Density scale.
	double P_SCALE; //!< Pressure scale.
	double E_SCALE; //!< Energy scale.
	static bool debugging; //!< Switch for outputting debug information to console.
	static bool freqPrinting; //!< Switch to printing data.
	static int percent; //!< Percent complete.
	static std::string progressMessage; //!< Prefix for progress bar.

	InputOutput(const PrintParameters& rp, const Scalings& sc);

	static void initProgressBar(std::string msg, MPIHandler& mpihandler);
	static void progressBar(const int& prcnow, const int& every, MPIHandler& mpihandler);
	static void endProgressBar(MPIHandler& mpihandler);

	//Output.
	void freqPrint(const Radiation* rad, Grid3D* gptr, MPIHandler& mpihandler);
	void printSTARBENCH(const int& i, const Radiation* rad, Grid3D* gptr, MPIHandler& mpih);
	void printBinary2D(const int& step, const double& t, const double& dt, Grid3D* gptr, MPIHandler& mpih) const;
	void print2D(const int& step, const double& t, const double& dt, Grid3D* gptr, MPIHandler& mpih) const;
	void printIF(const double& t, const Radiation* rad, Grid3D* gptr, MPIHandler& mpih) const;
	void printWeights(Grid3D* gptr);
	void printCellPathLength(Grid3D* gptr);
	//void printParams(double runtime, const GridParameters& gpar, const RadiationParameters& rpar,
		//const HydroParameters& hpar, const PrintParameters& ppar, const Scalings& scale) const;

	//Input.
	void fileToMap(const std::string& myString, std::map<double,double> myMap) const;

	//Formatting.
	std::string fortranformat(double value, int w, int d, int e);

	//Simulation timestep affecting methods.
	void addPrintTime(const double& t) {
		printTimes.push_back(t);
	}
	void reduceToPrint(const double& currTime, double& dt) {
		for (int i = 0; i < (int)printTimes.size(); i++) {
			if(currTime < printTimes[i] && currTime+dt > printTimes[i]) {
				dt = std::min(dt, printTimes[i]-currTime);
			}
		}
	}
};


#endif
