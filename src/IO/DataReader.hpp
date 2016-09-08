/** Provides the DataReader class.
 *
 * @file DataReader.hpp
 *
 * @author Harrison Steggles
 */

#ifndef DATAREADER_HPP_
#define DATAREADER_HPP_

#include <string>
#include <array>

class Fluid;
class DataParameters;

/**
 * @class DataPrinter
 *
 * @brief Handler for all data input operations.
 */
class DataReader {
public:
	static DataParameters readDataParameters(const std::string& filename);
	static void readGrid(const std::string& filename, const DataParameters& dp,  Fluid& fluid);
	static void patchGrid(const std::string& filename, const std::array<int, 3>& offset, Fluid& fluid );

private:

};



#endif /* DATAREADER_HPP_ */
