#include <iostream>
#include <map>

#include "parameterparser.h"
#include "parameters.h"
#include "tinyxml.h"

bool ParameterParser::parseParameters(const char* paramFile, IntegrationParameters& integParams,
		GridParameters& gridParams, HydroParameters& hydroParams, RadiationParameters& radParams, PrintParameters& printParams) {
	// Create XML document.
	TiXmlDocument xmlDoc;
	// Load the state file.
	if (!xmlDoc.LoadFile(paramFile)) {
		std::cerr << xmlDoc.ErrorDesc() << "\n";
		return false;
	}
	// Get the root element.
	TiXmlElement* pRoot = xmlDoc.RootElement();
	// Get the GRID root node and assign it to pGridRoot
	for (TiXmlElement* e = pRoot->FirstChildElement(); e != NULL; e = e->NextSiblingElement()) {
		if (e->Value() == std::string("INTEGRATION"))
			parseIntegrationParameters(e, integParams);
		if (e->Value() == std::string("GRID"))
			parseGridParameters(e, gridParams);
		if (e->Value() == std::string("RADIATION")) {
			parseRadiationParameters(e, radParams);
		}
		if (e->Value() == std::string("STARS")) {
			parseStarParameters(e, radParams);
		}
		if (e->Value() == std::string("HYDRODYNAMICS"))
			parseHydroParameters(e, hydroParams);
		if (e->Value() == std::string("PRINTING"))
			parsePrintParameters(e, printParams);
	}
	return true;
}

void ParameterParser::parseIntegrationParameters(TiXmlElement* pElement, IntegrationParameters& integParams) {
	TiXmlElement* pChild = pElement->FirstChildElement();
	pChild->Attribute("spatial_order", &integParams.ORDER_S);
	pChild->Attribute("temporal_order", &integParams.ORDER_T);
	pChild->Attribute("max_timestep", &integParams.DT_MAX);
}

void ParameterParser::parseGridParameters(TiXmlElement* pElement, GridParameters& gridParams) {
	TiXmlElement* pChild = pElement->FirstChildElement();
	pChild->Attribute("no_dimensions", &gridParams.ND);
	pChild->Attribute("no_cells_x", &gridParams.NCELLS[0]);
	pChild->Attribute("no_cells_y", &gridParams.NCELLS[1]);
	pChild->Attribute("no_cells_z", &gridParams.NCELLS[2]);
	pChild->Attribute("side_length", &gridParams.SIDE_LENGTH);

	std::map<std::string, Geometry> mGeometry;
	mGeometry["cartesian"] = CARTESIAN;
	mGeometry["cylindrical"] = CYLINDRICAL;
	mGeometry["spherical"] = SPHERICAL;
	if (mGeometry.find(pChild->Attribute("geometry")) != mGeometry.end())
		gridParams.GEOMETRY = mGeometry[pChild->Attribute("geometry")];

	std::map<std::string, Condition> mBC;
	mBC["free"] = FREE;
	mBC["reflecting"] = REFLECTING;
	mBC["outflow"] = OUTFLOW;
	mBC["inflow"] = INFLOW;
	if (mBC.find(pChild->Attribute("left_boundary_condition_x")) != mBC.end())
		gridParams.LBCondition[0] = mBC[pChild->Attribute("left_boundary_condition_x")];
	if (mBC.find(pChild->Attribute("left_boundary_condition_y")) != mBC.end())
		gridParams.LBCondition[1] = mBC[pChild->Attribute("left_boundary_condition_y")];
	if (mBC.find(pChild->Attribute("left_boundary_condition_z")) != mBC.end())
		gridParams.LBCondition[2] = mBC[pChild->Attribute("left_boundary_condition_z")];
	if (mBC.find(pChild->Attribute("right_boundary_condition_x")) != mBC.end())
		gridParams.RBCondition[0] = mBC[pChild->Attribute("right_boundary_condition_x")];
	if (mBC.find(pChild->Attribute("right_boundary_condition_y")) != mBC.end())
		gridParams.RBCondition[1] = mBC[pChild->Attribute("right_boundary_condition_y")];
	if (mBC.find(pChild->Attribute("right_boundary_condition_z")) != mBC.end())
		gridParams.RBCondition[2] = mBC[pChild->Attribute("right_boundary_condition_z")];
}

void ParameterParser::parseRadiationParameters(TiXmlElement* pElement, RadiationParameters& radParams) {
	TiXmlElement* pChild = pElement->FirstChildElement();
	pChild->Attribute("K1", &radParams.K1);
	pChild->Attribute("K2", &radParams.K2);
	pChild->Attribute("K3", &radParams.K3);
	pChild->Attribute("K4", &radParams.K4);
	pChild->Attribute("photoion_cross_section", &radParams.P_I_CROSS_SECTION);
	pChild->Attribute("case_b_recombination_coeff", &radParams.ALPHA_B);
	pChild->Attribute("tau_0", &radParams.TAU_0);
	pChild->Attribute("initial_no_density_hi", &radParams.NHI);
	pChild->Attribute("temperature_hi", &radParams.THI);
	pChild->Attribute("temperature_hii", &radParams.THII);

	std::string scheme = pChild->Attribute("integration_scheme");
	if (scheme == "implicit")
		radParams.SCHEME = IMPLICIT;
	else if (scheme == "explicit")
		radParams.SCHEME = EXPLICIT;
}

void ParameterParser::parseStarParameters(TiXmlElement* pElement, RadiationParameters& radParams) {
	for (TiXmlElement* pChild = pElement->FirstChildElement(); pChild != NULL; pChild = pChild->NextSiblingElement()) {
		StarParameters starParams;
		pChild->Attribute("cell_position_x", &starParams.POSITION[0]);
		pChild->Attribute("cell_position_y", &starParams.POSITION[1]);
		pChild->Attribute("cell_position_z", &starParams.POSITION[2]);
		parseBoolean(pChild->Attribute("snap_to_face_left_x"), starParams.FACE_SNAP[0]);
		parseBoolean(pChild->Attribute("snap_to_face_left_y"), starParams.FACE_SNAP[1]);
		parseBoolean(pChild->Attribute("snap_to_face_left_z"), starParams.FACE_SNAP[2]);
		parseBoolean(pChild->Attribute("snap_to_face_right_x"), starParams.FACE_SNAP[3]);
		parseBoolean(pChild->Attribute("snap_to_face_right_y"), starParams.FACE_SNAP[4]);
		parseBoolean(pChild->Attribute("snap_to_face_right_z"), starParams.FACE_SNAP[5]);
		pChild->Attribute("photon_energy", &starParams.PHOTON_ENERGY);
		pChild->Attribute("photon_rate", &starParams.PHOTON_RATE);
		radParams.vStarParams.push_back(starParams);
	}
}

void ParameterParser::parseHydroParameters(TiXmlElement* pElement, HydroParameters& hydroParams) {
	TiXmlElement* pChild = pElement->FirstChildElement();
	pChild->Attribute("gamma", &hydroParams.GAMMA);
	pChild->Attribute("density_floor", &hydroParams.DFLOOR);
	pChild->Attribute("pressure_floor", &hydroParams.PFLOOR);
}

void ParameterParser::parsePrintParameters(TiXmlElement* pElement, PrintParameters& printParams) {
	TiXmlElement* pChild = pElement->FirstChildElement();
	parseBoolean(pChild->Attribute("print_solution"), printParams.PRINT2D_ON);
	parseBoolean(pChild->Attribute("print_ionisation_front"), printParams.PRINTIF_ON);
}

void ParameterParser::parseBoolean(std::string str, bool& variable) {
	if (str == std::string("true"))
		variable = true;
	else if (str == std::string("false"))
		variable = false;
}


