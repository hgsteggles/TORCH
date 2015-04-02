/** Provides the DataPrinter class.
 *
 * @file DataPrinter.hpp
 *
 * @author Harrison Steggles
 *
 * @date 13/01/2014 - The first version.
 * @date 04/01/2014 - Arguments now passed by const reference when appropriate.
 * @date 12/02/2014 - Progress bar added which outputs progress to command line and sets the frequency for
 * printing data, initProgressBar needs to be called first and progressBar() needs to be called at the end.
 * @date 21/07/2014 - Replaced const reference to int/double with a const copy.
 * @date 24/08/2014 - Renamed to DataPrinter.
 * @date 03/12/2014 - No longer uses boost::iostreams. Created zlib wrapper class (to reduce dependencies).
 * @date 03/12/2014 - Throws runtime_errors instead of outputting to console.
 */

#ifndef DATAPRINTER_HPP_
#define DATAPRINTER_HPP_

#include <map>
#include <memory>
#include <string>
#include <vector>

class PrintParameters;
class Converter;
class Fluid;
class Radiation;
class Hydrodynamics;
class Constants;
class Fluid;
class Grid;

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
	void printMinMax(const std::string& filename, const Grid& grid) const;
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

#endif // DATAPRINTER_HPP_
