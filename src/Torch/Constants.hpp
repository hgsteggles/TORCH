/** Provides Constants and EnumParser classes.
 * @file Constants.hpp
 *
 *  @author "Harrison Steggles"
 *
 *  @date 24/11/2014 - the first version.
 */

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include "Converter.hpp"
#include "Common.hpp"

#include <string>
#include <map>
#include <stdexcept>

class Constants;

/**
 * @class EnumParser
 *
 * @brief Parses strings and converts them into enums.
 *
 * @version 0.8, 24/11/2014
 */
template <typename T> class EnumParser {
	friend class Constants;
public:
	EnumParser(){};

	T parseEnum(const std::string &value) {
		typename std::map<std::string, T>::const_iterator iValue = enumMap.find(value);
		if (iValue  == enumMap.end()) {
			typename std::map<std::string, T>::const_iterator dValue = enumMap.find("default");
			if (dValue == enumMap.end())
				throw std::runtime_error("EnumParser: invalid key [" + value + "]");
			return dValue->second;
		}
		return iValue->second;
	}
private:
	std::map<std::string, T> enumMap;
};

/**
 * @class Constants
 *
 * @brief Contains all constants that may be needed at any place in the code. Also has a Converter handy for unit conversions.
 *
 * @version 0.8, 24/11/2014
 */
class Constants {
public:
	Constants();
	void initialise_MLT(double mass, double length, double time);
	void initialise_DPT(double density, double pressure, double time);

	Converter converter;

	double specificGasConstant = 0; //!< Specific Gas Constant.
	double rydbergEnergy = 0; //!< Rydberg Energy.
	double boltzmannConst = 0; //!< Boltzmann's Constant.
	double dustExtinctionCrossSection = 0; //!< Dust Extinction Cross Section (Baldwin et. al. 1991) [cm2 H-1].
	double hydrogenMass = 0; //!< Mass of hydrogen in g.
	double pi = 0; //!< PI.

	int nd = 3; //!< Number of Dimensions.
	double dfloor = 0;
	double pfloor = 0;
	double tfloor = 0;

	//Voronov (1997,ADANDT,65,1).
	double voronov_A = 0; //!< Constant A from Voronov (1997,ADANDT,65,1).
	double voronov_X = 0; //!< Constant X from Voronov (1997,ADANDT,65,1).
	double voronov_K = 0; //!< Constant K from Voronov (1997,ADANDT,65,1).
	double voronov_P = 0; //!< Constant P from Voronov (1997,ADANDT,65,1).
	double voronov_DE = 0; //!< Constant DE from Voronov (1997,ADANDT,65,1).

	EnumParser<Geometry> geometryParser;
	EnumParser<Condition> conditionParser;
	EnumParser<Scheme> schemeParser;
	EnumParser<Coupling> couplingParser;


private:
	void initialise();
};

#endif // CONSTANTS_HPP_
