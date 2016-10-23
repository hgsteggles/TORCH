/** Provides Grid typedefs and enums.
 *
 * @file Common.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef COMMON_HPP_
#define COMMON_HPP_

#include <stddef.h>
#include <array>

class Boundary;
class GridCell;
class GridJoin;

template <class T, size_t ROW, size_t COL>
using Array2D = std::array<std::array<T, COL>, ROW>;

struct UID {
	enum ID {DEN, PRE, HII, ADV, VEL, N=7};
};
struct RID {
	enum ID {HII_A, TAU, TAU_A, DTAU, DTAU_A, HEAT, N};
};
struct TID {
	enum ID {COL_DEN, DCOL_DEN, HEAT, RATE, N};
};
struct HID {
	enum ID {IMLC, NMLC, RHII, CEHI, CIEC, NMC, EUVH, FUVH, IRH, CRH, TOT, N};
};

enum class Geometry : unsigned int {CARTESIAN, CYLINDRICAL, SPHERICAL};
enum class Condition : unsigned int {FREE, REFLECTING, OUTFLOW, INFLOW, PERIODIC, PARTITION};
enum class Scheme : unsigned int {IMPLICIT, IMPLICIT2, EXPLICIT};
enum class Coupling : unsigned int {TWO_TEMP_ISOTHERMAL, NON_EQUILIBRIUM, OFF};

using FluidArray = std::array<double, UID::N>;
using RadArray = std::array<double, RID::N>;
using ThermoArray = std::array<double, TID::N>;
using HeatArray = std::array<double, HID::N>;
using Vec3 = std::array<double, 3>;
using Coords = std::array<int, 3>;

using Vec3 = std::array<double, 3>;
using Coords = std::array<int, 3>;
using ConditionVec = std::array<Condition, 3>;

#endif // COMMON_HPP_
