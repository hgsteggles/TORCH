/** Provides the DataPrinter class.
 *
 * @file DataPrinter.hpp
 *
 * @author Harrison Steggles
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
 */
class DataPrinter {
public:
	void initialise(std::shared_ptr<Constants> c, std::string output_directory);

	//Output.
	void freqPrint(const Radiation& rad, const Grid& grid) const;
	void printSTARBENCH(const Radiation& rad, const Hydrodynamics& hydro, Fluid& fluid);
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
