/** Provides the Converter class.
 * @file Converter.hpp
 *
 * @author Harrison Steggles
 *
 * @date 24/11/2014 - the first version.
 */

#ifndef CONVERTER_HPP_
#define CONVERTER_HPP_

/**
 * @class Converter
 *
 * @brief Converts values to/from code units from/to cgs real units.
 *
 * @version 0.8, 24/11/2014
 */
class Converter {
public:
	Converter();
	void set_mass_length_time(const double mass, const double length, const double time);
	void set_rho_pressure_time(const double rho, const double pressure, const double time);
	double toCodeUnits(const double val, const double mass_index, const double length_index, const double time_index) const;
	double fromCodeUnits(const double val, const double mass_index, const double length_index, const double time_index) const;

	double EV_2_ERGS(double val_in_ev) const;
	double ERG_2_EV(double val_in_yr) const;
	double YR_2_S(double val_in_yr) const;
	double S_2_YR(double val_in_yr) const;
	double PC_2_CM(double val_in_yr) const;
	double CM_2_PC(double val_in_yr) const;
	double MO_2_G(double val_in_yr) const;
private:
	double M = 0, L = 0, T = 0;
	double convertCodeUnits(const double val, const double mass_index, const double length_index, const double time_index, const bool& from) const;
};

#endif // CONVERTER_HPP_
