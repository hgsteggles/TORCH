/**
 * Provides global constants and enums.
 * @file constants.hpp
 *
 *  @author "Harrison Steggles"
 *  @date 05/02/2014 - Constants and Enums all in their own file.
 */

#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

enum Geometry {CARTESIAN, CYLINDRICAL, SPHERICAL};
enum Condition {FREE, REFLECTING, OUTFLOW, INFLOW};
enum Scheme {IMPLICIT, IMPLICIT2, EXPLICIT};
enum SendID {PARTITION_MSG, RADIATION_MSG, PRINT2D_MSG, CFL_COLLECT, CFL_BROADCAST, PRINTIF_NEXT_MSG, PRINTIF_FOUND_MSG, PRINTIF_IF_MSG, PRINTSTARBENCH_MSG};
const int NU = 6;
const int iden = 0;
const int ivel = 1;
const int ipre = 4;
const int ihii = 5;
const int NR = 7;
const int ihiita = 0;
const int itau = 1;
const int itauta = 2;
const int idtau = 3;
const int idtauta = 4;
const int icolden = 5;
const int icool = 6;

extern bool SWEEPX;
/* MICROPHYSICS */
const bool RTCOUPLING = true;
const bool MONOCHROME = true;
const bool RECOMBINATIONS = true;
const bool COLLISIONS = false;
const bool COOLING = false;
/* CONVERSIONS */
const double YR2S = 3.15569e7;
const double S2YR = 1.0/YR2S;
const double PC2CM = 3.09e18;
const double CM2PC = 1.0/PC2CM;
const double MO2G = 2e33;
const double EV2ERG = 1.60217646e-12;
const double ERG2EV = 1.0/EV2ERG;
/* CONSTANTS */
const double PI = 3.14159265359;
const double GAS_CONST = 8.314462e7; // ergs mol^-1 K*-1
const double BOLTZMANN = 1.3806488e-16; //erg K^-1
const double H_MASS_G = 1.674e-24; // g


#endif /* CONSTANTS_HPP_ */
