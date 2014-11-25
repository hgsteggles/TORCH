#include "Constants.hpp"

Constants::Constants() {
	geometryParser.enumMap["cartesian"] = Geometry::CARTESIAN;
	geometryParser.enumMap["cylindrical"] = Geometry::CYLINDRICAL;
	geometryParser.enumMap["spherical"] = Geometry::SPHERICAL;
	geometryParser.enumMap["default"] = Geometry::CARTESIAN;

	conditionParser.enumMap["free"] = Condition::FREE;
	conditionParser.enumMap["reflecting"] = Condition::REFLECTING;
	conditionParser.enumMap["outflow"] = Condition::OUTFLOW;
	conditionParser.enumMap["inflow"] = Condition::INFLOW;
	conditionParser.enumMap["periodic"] = Condition::PERIODIC;
	conditionParser.enumMap["default"] = Condition::FREE;

	schemeParser.enumMap["implicit"] = Scheme::IMPLICIT;
	schemeParser.enumMap["implicit2"] = Scheme::IMPLICIT2;
	schemeParser.enumMap["explicit"] = Scheme::EXPLICIT;
	schemeParser.enumMap["default"] = Scheme::IMPLICIT;

	couplingParser.enumMap["tti"] = Coupling::TWO_TEMP_ISOTHERMAL;
	couplingParser.enumMap["neq"] = Coupling::NON_EQUILIBRIUM;
	couplingParser.enumMap["off"] = Coupling::OFF;
	couplingParser.enumMap["default"] = Coupling::OFF;
}

void Constants::initialise() {
	hydrogenMass = converter.toCodeUnits(1.674e-24, 1, 0, 0);
	specificGasConstant = converter.toCodeUnits(8.314462e7, 0, 2, -2);
	boltzmannConst = converter.toCodeUnits(1.3806488e-16, 1, 2, -2);
	rydbergEnergy = converter.toCodeUnits(converter.EV_2_ERGS(13.6), 1, 2, -2);
	dustExtinctionCrossSection = converter.toCodeUnits(5.0e-22, 0, 2, 0);
	pi = 3.14159265359;

	voronov_A = converter.toCodeUnits(2.91e-8, 0, 3, -1);
	voronov_X = 0.232;
	voronov_K = 0.39;
	voronov_P = 0.0;
	voronov_DE = rydbergEnergy/boltzmannConst;
}

void Constants::initialise_MLT(double mass, double length, double time) {
	converter.set_mass_length_time(mass, length, time);
	initialise();
}

void Constants::initialise_DPT(double density, double pressure, double time) {
	converter.set_rho_pressure_time(density, pressure, time);
	initialise();
}


