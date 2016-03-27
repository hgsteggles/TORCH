/*
 * DataReader.hpp
 *
 *  Created on: 4 Mar 2015
 *      Author: harry
 */

#ifndef DATAREADER_HPP_
#define DATAREADER_HPP_

#include <string>
#include <array>

class Fluid;
class DataParameters;

class DataReader {
public:
	static DataParameters readDataParameters(const std::string& filename);
	static void readGrid(const std::string& filename, const DataParameters& dp,  Fluid& fluid);
	static void patchGrid(const std::string& filename, const std::array<int, 3>& offset, Fluid& fluid );

private:

};



#endif /* DATAREADER_HPP_ */
