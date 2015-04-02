#include "Converter.hpp"

#include <cmath>
#include <iostream>


Converter::Converter() : M(MO_2_G(1.0)), L(PC_2_CM(14.6)), T(YR_2_S(10000.0)) {

}

void Converter::set_mass_length_time(const double mass, const double length, const double time) {
	M = (mass > 0) ? mass : M;
	L = (length > 0) ? length : L;
	T = (time > 0) ? time : T;
}

void Converter::set_rho_pressure_time(const double rho, const double pressure, const double time) {
	double rho_0 = (rho > 0) ? rho : M/(L*L*L);
	double pressure_0 = (pressure > 0) ? pressure : rho_0*L*L/(T*T);
	double time_0 = (time > 0) ? time : T;

	T = time_0;
	M = std::sqrt((pressure_0/rho_0)*pressure_0*pressure_0)*T*T*T;
	L = std::sqrt(pressure_0/rho_0)*T;
}

double Converter::toCodeUnits(const double val, const double mass_index, const double length_index, const double time_index) const {
	return convertCodeUnits(val, mass_index, length_index, time_index, true);
}

double Converter::fromCodeUnits(const double val, const double mass_index, const double length_index, const double time_index) const {
	return convertCodeUnits(val, mass_index, length_index, time_index, false);
}

double Converter::convertCodeUnits(const double val, const double mass_index, const double length_index, const double time_index,
		const bool& to_code_units) const {
	double MM, LL, TT;
	MM = mass_index >= 0 ? M : 1.0/M;
	MM = to_code_units ? 1.0/MM : MM;
	LL = length_index >= 0 ? L : 1.0/L;
	LL = to_code_units ? 1.0/LL : LL;
	TT = time_index >= 0 ? T : 1.0/T;
	TT = to_code_units ? 1.0/TT : TT;
	int im = (int)std::abs(mass_index);
	int il = (int)std::abs(length_index);
	int it = (int)std::abs(time_index);
	double result = val;
	while (im >= 1 || il >= 1 || it >= 1) {
		if (im-- >= 1) result *= MM;
		if (il-- >= 1) result *= LL;
		if (it-- >= 1) result *= TT;
	}
	double left = std::abs(mass_index) - (int)std::abs(mass_index);
	if (left != 0)
		result *= std::pow(MM, left);
	left = std::abs(length_index) - (int)std::abs(length_index);
	if (left != 0)
		result *= std::pow(LL, left);
	left = std::abs(time_index) - (int)std::abs(time_index);
	if (left != 0)
		result *= std::pow(TT, left);
	return result;
}

double Converter::EV_2_ERGS(double val_in_ev) const {
	return val_in_ev*1.60217646e-12;
}

double Converter::ERG_2_EV(double val_in_erg)  const{
	return val_in_erg/EV_2_ERGS(1.0);
}

double Converter::YR_2_S(double val_in_yr) const {
	return val_in_yr*3.15569e7;
}

double Converter::S_2_YR(double val_in_s) const {
	return val_in_s/YR_2_S(1.0);
}

double Converter::PC_2_CM(double val_in_pc) const {
	return val_in_pc*3.09e18;
}

double Converter::CM_2_PC(double val_in_cm) const {
	return val_in_cm/PC_2_CM(1.0);
}

double Converter::MO_2_G(double val_in_mo) const {
	return val_in_mo*2.0e33;
}

void Converter::printInfo() {
	std::cout << "M = " << M << "\n";
	std::cout << "L = " << L << "\n";
	std::cout << "T = " << T << "\n";
}
