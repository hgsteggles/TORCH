/** Provides the DataPrinter class.
 *
 * @file DataPrinter.hpp
 *
 * @author Harrison Steggles
 *
 * @date 13/01/2014 - the first version.
 * @date 04/01/2014 - arguments now passed by const reference when appropriate.
 * @date 12/02/2014 - progress bar added which outputs progress to command line and sets the frequency for
 * printing data, initProgressBar needs to be called first and progressBar() needs to be called at the end.
 * @date 21/07/2014 - replaced const reference to int/double with a const copy.
 * @date 24/08/2014 - renamed to DataPrinter.
 */

#ifndef DATAPRINTER_HPP_
#define DATAPRINTER_HPP_

#include "MPI_Wrapper.hpp"

#include <map>
#include <vector>
#include <memory>

class PrintParameters;
class Converter;
class Grid;
class Radiation;
class Hydrodynamics;
class Constants;
class Fluid;

/**
 * @class DataPrinter
 *
 * @brief Handler for all data output operations.
 *
 * @version 0.8, 24/11/2014
 */
class DataPrinter {
public:
	void initialise(std::shared_ptr<Constants> c, std::string output_directory);

	//Output.
	void freqPrint(const Radiation& rad, const Grid& grid) const;
	void printSTARBENCH(const Radiation& rad, const Hydrodynamics& hydro, const Fluid& fluid);
	void printBinary2D(const int step, const double t, const Grid& grid) const;
	void print2D(const std::string& append_name, const double t, const Grid& grid) const;
	void printHeating(const int step, const double t, const Grid& grid) const;
	void printVariables(const int step, const double t, const Grid& grid) const;
	void printVariable(const int step, const double t, const Grid& grid) const;
	void printWeights(const Grid& grid) const;

	//Input.
	void fileToMap(const std::string& myString, std::map<double,double> myMap) const;

	//Formatting.
	std::string fortranformat(double value, int w, int d, int e) const;

	//Simulation timestep affecting methods.
	void addPrintTime(const double t);
	void reduceToPrint(const double currTime, double& dt) const;

private:
	std::shared_ptr<Constants> consts = nullptr;
	std::string dir2D = "tmp/";
	std::vector<double> printTimes;
	std::vector<bool> printDone;
	bool printing_on = true;
};

struct PrintParameters {
	PrintParameters();
	std::string dir2D = "tmp/";
	void printInfo() const;
};

#endif /* DATAPRINTER_HPP_ */
