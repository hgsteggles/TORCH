/**
 * 
 * @file parameterparser.h
 *
 *  @author "Harrison Steggles"
 *  @date
 */

#ifndef PARAMETERPARSER_HPP_
#define PARAMETERPARSER_HPP_

#include <iostream>
#include <vector>
#include "tinyxml.h"

class IntegrationParameters;
class GridParameters;
class HydroParameters;
class RadiationParameters;
class PrintParameters;

class ParameterParser {
public:
	bool parseParameters(const char* stateFile, IntegrationParameters& integParams, GridParameters& gridParams,
			HydroParameters& hydroParams, RadiationParameters& radParams, PrintParameters& printParams);
private:
	void parseIntegrationParameters(TiXmlElement* pElement, IntegrationParameters& integParams);
	void parseGridParameters(TiXmlElement* pElement, GridParameters& gridParams);
	void parseHydroParameters(TiXmlElement* pElement, HydroParameters& hydroParams);
	void parseRadiationParameters(TiXmlElement* pElement, RadiationParameters& radParams);
	void parseStarParameters(TiXmlElement* pElement, RadiationParameters& radParams);
	void parsePrintParameters(TiXmlElement* pElement, PrintParameters& printParams);
	void parseBoolean(std::string str, bool& variable);
};

#endif /* PARAMETERPARSER_HPP_ */
