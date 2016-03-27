/** Provides the parameter classes.
 *
 * @file Parameters.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef PARAMETERS_HPP_
#define PARAMETERS_HPP_

#include <array>
#include <memory>
#include <string>

class Constants;

struct FluidParameters;
struct GridParameters;
struct RadiationParameters;
struct ThermoParameters;
struct StarParameters;

struct TorchParameters {
	std::string setupFile = "";
	std::string initialConditions = "";
	int nd = 0; //!< Number of dimensions.
	double sideLength = 0; //!< The side length of the simulation line/square/cube.
	std::array<int, 3> ncells = std::array<int, 3>{ 0, 0, 0 }; //!< Array holding the number of grid cells along each dimension.
	std::array<int, 3> coreCells = std::array<int, 3>{ 0, 0, 0 }; //!< Array holding the number of grid cells along in each dimension in each processor.
	std::array<std::string, 3> leftBC = std::array<std::string, 3>{ "free", "free", "free" }; //!< Array of left boundary conditions for each dimension.
	std::array<std::string, 3> rightBC = std::array<std::string, 3>{ "free", "free", "free" }; //!< Array of right boundary conditions for each dimension.
	std::string geometry = "cartesian"; //!< The GEOMETRY of the grid [CARTESIAN, CYLINDRICAL, SPHERICAL].

	std::string patchfilename = "";
	std::array<int, 3> patchoffset = std::array<int, 3>{ 0, 0, 0 };

	std::string riemannSolver = "hll";
	std::string slopeLimiter = "falle";
	std::string rt_scheme = "implicit";  //!< Ionisation fraction integration scheme.
	std::string rt_coupling = "off";

	int spatialOrder = 0;
	int temporalOrder = 0;
	double tmax = 0;
	double dt_max = 0;
	bool radiation_on = false;
	bool cooling_on = false;
	bool debug = true;
	std::string outputDirectory = "tmp/";

	double dfloor = 0;
	double pfloor = 0;
	double tfloor = 0;
	double dscale = 0;
	double pscale = 0;
	double tscale = 0;

	double K1 = 0;
	double K2 = 0;
	double K3 = 0;
	double K4 = 0; //!< [Mackey 2012 (table 1)] timestep constant.
	double photoIonCrossSection = 0; //!< Photoionisation cross-section (cm2.scale-1)
	double alphaB = 0; //!< Case B radiative recombination rate coefficient (cm3.t-1.scale-1).
	double tau0 = 0; //!< Min. tau in nearest neighbour weights [Mellema et. al. 2006 (A.5)].
	double minX = 0; //!< Minimum fraction of ionised hydrogen.
	double THI = 0; //!< Temperature fix for fully neutral gas.
	double THII = 0; //!< Temperature fix for fully ionized gas.
	bool collisions_on = false; //!< Include collisional ionizations.

	bool star_on = false;
	std::array<int, 3> star_position = std::array<int, 3>{ 0, 0, 0 };
	std::array<bool, 3> faceSnap = std::array<bool, 3>{ false, false, false };
	double photonEnergy = 0;
	double photonRate = 0;
	double massLossRate = 0;
	double windVelocity = 0;
	double windTemperature = 0;
	int windCellRadius = 0;

	double thermoHII_Switch = 0;
	double heatingAmplification = 1.0; //!< Heating amplification/reduction hack.
	bool thermoSubcycling = true;
	bool minTempInitialState = false;

	double massFractionH = 0; //!< Mass fraction of hydrogen.
	double heatCapacityRatio = 0;

	void initialise(std::shared_ptr<Constants>& consts);
	GridParameters getGridParameters();
	FluidParameters getFluidParameters();
	RadiationParameters getRadiationParameters();
	ThermoParameters getThermoParameters();
	StarParameters getStarParameters();
};

struct FluidParameters {
	double heatCapacityRatio = 0;
	double massFractionH = 0;
};

struct GridParameters {
	std::array<int, 3> ncells; //!< Array holding the number of grid cells along each dimension.
	int spatialOrder;
	double sideLength; //!< The side length of the simulation line/square/cube.
	std::array<std::string, 3> leftBC; //!< Array of left boundary conditions for each dimension.
	std::array<std::string, 3> rightBC; //!< Array of right boundary conditions for each dimension.
	std::string geometry; //!< The GEOMETRY of the grid [CARTESIAN, CYLINDRICAL, SPHERICAL].
};

struct DataParameters {
	double time = 0;
	std::array<int, 3> ncells = std::array<int, 3>{ 0, 0, 0 }; //!< Array holding the number of grid cells along each dimension.
	int nd = 0;
	double dx = 0;
	double sideLength = 0; //!< The side length of the simulation line/square/cube.
};

struct RadiationParameters {
	double K1 = 0; //!< [Mackey 2012 (table 1)] timestep constant.
	double K2 = 0; //!< [Mackey 2012 (table 1)] timestep constant.
	double K3 = 0; //!< [Mackey 2012 (table 1)] timestep constant.
	double K4 = 0; //!< [Mackey 2012 (table 1)] timestep constant.
	double photoIonCrossSection = 0; //!< Photoionisation cross-section (cm2.scale-1)
	double alphaB = 0; //!< Case B radiative recombination rate coefficient (cm3.t-1.scale-1).
	double tau0 = 0; //!< Min. tau in nearest neighbour weights [Mellema et. al. 2006 (A.5)].
	double minX = 0;
	double THI = 0; //!< Temperature fix for fully neutral gas.
	double THII = 0; //!< Temperature fix for fully ionized gas.
	std::string scheme = "implicit";  //!< Ionization fraction integration scheme.
	double massFractionH = 0; //!< Mass fraction of hydrogen.
	double heatingAmplification = 0;
	bool collisions_on = false; //!< Include collisional ionizations.
	std::string coupling = "off";
};

struct ThermoParameters {
	double thermoHII_Switch = 0;
	double heatingAmplification = 1.0; //!< Heating amplification/reduction hack.
	double massFractionH = 0;
	bool thermoSubcycling = true;
	bool minTempInitialState = false;
};

struct StarParameters {
	bool on = false;
	std::array<int, 3> position = std::array<int, 3>{ 0, 0, 0 };
	std::array<bool, 3> faceSnap = std::array<bool, 3>{ false, false, false };
	double photonEnergy = 0;
	double photonRate = 0;
	int windCellRadius = 0;
	double massLossRate = 0;
	double windVelocity = 0;
	double windTemperature = 0;
};


#endif // PARAMETERS_H_
