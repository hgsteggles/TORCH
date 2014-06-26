#include "thermodynamics.h"
#include "grid3d.h"
#include "gridcell.h"
#include "radiation.h"
#include "star.h"
#include "parameters.h"
#include "recipes.h"
#include "constants.h"

Thermodynamics::Thermodynamics(const Scalings& scale, Grid3D& g3d) : grid(g3d) {
	spec_gas_const = scale.toCodeUnits(8.314462e7, 0, 2, -2); // g-1 ergs mol^-1 K*-1
	kB = scale.toCodeUnits(1.3806488e-16, 1, 2, -2);
	sigmaV = 5.0e-22; //(Baldwin et al. 1991)
	z0 = 5.0e-4;
	excessEnergy = scale.toCodeUnits(5*EV2ERG, 1, 2, -2);
	T1 = 33610;
	T2 = 2180;
	T3 = 28390;
	T4 = 1780;
	imlc = scale.toCodeUnits(2.905e-19, 1.0, 5.0, -3.0);
	nmlc = scale.toCodeUnits(4.477e-20, 1.0, 5.0, -3.0);
	ciec_minT = 5.0e4;
	ciec = scale.toCodeUnits(3.485e-15, 1.0, 5.0, -3.0);
	cxhi_damp = 5.0e5;
	n0 = scale.toCodeUnits(1.0e6, 1.0, 2.0, -2.0);
	nmc = scale.toCodeUnits(3.981e-27, 1.0, 3.8, -3.0);
	//fuvh_a = scale.toCodeUnits(1.9e-26, 1.0, 5.0, -2.0);
	//fuvh_b = scale.toCodeUnits(1.00000, 0.0, 1.0, -1.0);
	//fuvh_c = 6.4;
	fuvh_a = scale.toCodeUnits(1.9e-26, 1.0, 4.0, -2.0);
	fuvh_b = scale.toCodeUnits(1.00000, 0.0, 0.0, -1.0);
	fuvh_c = scale.toCodeUnits(6.40000, 0.0, -1.0, 0.0);
	irh_a = scale.toCodeUnits(7.7e-32, 1.0, 4.0, -2.0);
	irh_b = scale.toCodeUnits(3.0e4, 0.0, -3.0, 0.0);
	crh = scale.toCodeUnits(5.0e-27, 1.0, 2.0, -3.0);

	/*
	std::cout << "imlc = " << imlc << "\n";
	std::cout << "nmlc = " << nmlc << "\n";
	std::cout << "ciec = " << ciec << "\n";
	std::cout << "n0 = " << n0 << "\n";
	std::cout << "nmc = " << nmc << "\n";
	std::cout << "fuvh_a = " << fuvh_a << "\n";
	std::cout << "fuvh_b = " << fuvh_b << "\n";
	std::cout << "fuvh_c = " << fuvh_c << "\n";
	std::cout << "irh_a = " << irh_a << "\n";
	std::cout << "irh_b = " << irh_b << "\n";
	std::cout << "crh = " << crh << "\n";
	exit(EXIT_FAILURE);
	*/

	initCollisionalExcitationHI(scale);
	initRecombinationHII(scale);
}

/**
 * @brief Cubic spline fit for the Aggarwal (1993) collisional excitation of HI data.
 * Tabulated values for the cooling rate from Raga, Mellema, \& Lundqvist, 1997, ApJS, 109, 517. The rates are for collisionally
 * excited cooling from neutral hydrogen (H I), and the rates are per electron, per HI atom, in [erg.cm3/s]. Spline is fit in
 * log-log space and the slopes off the end of the fit are also logarithmic.
 */
void Thermodynamics::initCollisionalExcitationHI(const Scalings& scale) {

	double T[26] = {3162.2776602, 3981.0717055, 5011.8723363, 6309.5734448, 7943.2823472,
			10000.0000000, 12589.2541179, 15848.9319246, 19952.6231497, 25118.8643151,
			31622.7766017, 39810.7170553, 50118.7233627, 63095.7344480, 79432.8234724,
			100000.0000000, 125892.5411794, 158489.3192461, 199526.2314969, 251188.6431510,
			316227.7660168, 398107.1705535, 501187.2336273, 630957.3444802, 794328.2347243,
			1000000.0000000};

	double R[26] = {1.150800e-34, 2.312065e-31, 9.571941e-29, 1.132400e-26,	4.954502e-25,
			9.794900e-24, 1.035142e-22, 6.652732e-22, 2.870781e-21, 9.036495e-21, 2.218196e-20,
			4.456562e-20, 7.655966e-20, 1.158777e-19, 1.588547e-19, 2.013724e-19, 2.393316e-19,
			2.710192e-19, 2.944422e-19, 3.104560e-19, 3.191538e-19, 3.213661e-19, 3.191538e-19,
			3.126079e-19, 3.033891e-19, 2.917427e-19};

	for (int i = 0; i < 26; i++) {
		cxhi_T.push_back(std::log10(T[i]));
		cxhi_rate.push_back(std::log10(scale.toCodeUnits(R[i], 1, 5, -3)));
	}
	cxhi_rate2 = NumericalRecipes::spline(cxhi_T, cxhi_rate, 1.e99, 1.e99);
	// Logarithmic slopes at either end of the domain.
	cxhi_minSlope = (cxhi_rate[1]-cxhi_rate[0])/(cxhi_T[1]-cxhi_T[0]);
	cxhi_maxSlope = (cxhi_rate[cxhi_T.size()-1]-cxhi_rate[cxhi_T.size()-2])/(cxhi_T[cxhi_T.size()-1]-cxhi_T[cxhi_T.size()-2]);
}

/**
 * @brief Cubic spline fit for Hummer (1994) HII recombination cooling rate data.
 */
void Thermodynamics::initRecombinationHII(const Scalings& scale) {
	double coolb[31] = {8.287e-11, 7.821e-11, 7.356e-11, 6.892e-11, 6.430e-11, 5.971e-11,
					5.515e-11, 5.062e-11, 4.614e-11, 4.170e-11, 3.734e-11, 3.306e-11, 2.888e-11,
					2.484e-11, 2.098e-11, 1.736e-11, 1.402e-11, 1.103e-11, 8.442e-12, 6.279e-12,
					4.539e-12, 3.192e-12, 2.185e-12, 1.458e-12, 9.484e-13, 6.023e-13, 3.738e-13,
					2.268e-13, 1.348e-13, 7.859e-14, 4.499e-14};
	for (int i = 0; i < 31; i++) {
		hr_T.push_back(std::exp(std::log(10.0)*(1.0 +0.2*static_cast<double>(i))));
		hr_rate.push_back(scale.toCodeUnits(coolb[i]/std::sqrt(hr_T[i]), 0, 3, -1));
	}
	hr_rate2 = NumericalRecipes::spline(hr_T, hr_rate, 1.e99,1.e99);
	hr_minSlope = (std::log10(hr_rate[1])-std::log10(hr_rate[0]))/(std::log10(hr_T[1])-std::log10(hr_T[0]));
	hr_maxSlope = (std::log10(hr_rate[30])-std::log10(hr_rate[29]))/(std::log10(hr_T[30])-std::log10(hr_T[29]));
}

double Thermodynamics::fluxFUV(const double& Q_FUV, const double dist_sqrd) {
	if (dist_sqrd != 0)
		return Q_FUV/(1.2e7*4*PI*dist_sqrd);
	else
		return 0;
}

/**
 * @brief Cubic spline interpolation of the collisional excitation cooling rate of HI.
 * The spline is fit in log-log space, and the slopes off the end of the fit are also logarithmic, so we take the log of T, get
 * log10(rate), and then return exp10() of the rate.
 * @param T Gas temperature.
 * @return Collisional excitation cooling rate of HI.
 */
double Thermodynamics::collisionalExcitationHI(const double& nH, const double& HIIFRAC, const double& T) {
	double rate = 0.0;
	double logT = log10(T);

	if (logT < cxhi_T[0])
		rate = cxhi_rate[0] + cxhi_minSlope*(logT - cxhi_T[0]);
	else if (logT > cxhi_T[cxhi_T.size()-1])
		rate = cxhi_rate[cxhi_T.size()-1] + cxhi_maxSlope*(logT - cxhi_T[cxhi_T.size()-1]);
	else
		rate = NumericalRecipes::splint(cxhi_T, cxhi_rate, cxhi_rate2, logT);

	return HIIFRAC*(1.0-HIIFRAC)*nH*nH*std::exp((2.302585093*rate)-((T/cxhi_damp)*(T/cxhi_damp)));
}

/**
 * Cubic spline interpolation of the recombination cooling rate of HII.
 * @param T Gas temperature.
 * @return Recombination cooling rate of HII.
 */
double Thermodynamics::recombinationHII(const double& nH, const double& HIIFRAC, const double& T) {
	double rate = 0;
	if (T > hr_T[hr_T.size()-1])
		rate = hr_rate[hr_T.size()-1] *std::pow(T/hr_T[hr_T.size()-1], hr_maxSlope);
	else if ( T < hr_T[0])
		rate = hr_rate[0]*std::pow(T/hr_T[0], hr_minSlope);
	else
		rate = NumericalRecipes::splint(hr_T, hr_rate, hr_rate2, T);

	return HIIFRAC*HIIFRAC*nH*nH*kB*T*rate;
}

/**
 * @brief Calculates the cooling rate of gas in a grid cell due to atomic processes.
 * Cooling due to collisionally excited optical lines of ionized metals; collisionally
 * excited lines of neutral metals; free-free and free-bound transitions of ionized
 * hydrogen; collisionally excited lines of neutral hydrogen; collisional ionization
 * equilibrium-cooling; and CLOUDY PDR models.
 * @param nH Hydrogen number density.
 * @param HIIFRAC Ionized hydrogen fraction.
 * @param T Gas temperature.
 * @return Cooling rate.
 */
double Thermodynamics::coolingProcesses(const double& nH, const double& HIIFRAC, const double& T) {
	double ne = nH*(HIIFRAC); //Ionized hydrogen no. density.
	double nn = nH*(1.0-HIIFRAC); //Neutral hydrogen no. density.
	double rate = 0;
	//Ionised metal line cooling (Henney et al. 2009, eq. A9).
	rate += imlc*z0*ne*ne*std::exp(-T1/T - (T2/T)*(T2/T));

	//Neutral metal line cooling (Henney et al. 2009, eq. A10)
	rate += nmlc*z0*ne*nn*std::exp(-T3/T - (T4/T)*(T4/T));
	//Free–free and free–bound transitions of ionized hydrogen (Henney et al. 2009, eq. A11).
	//rate += ne*ne*(1.34e-11*BOLTZMANN*sqrt(T)); //Osterbrock's (1989) free-free H+ function, Hummer says good to 30% in worst case.
	rate += recombinationHII(nH, HIIFRAC, T);
	//Collisionally excited lines of neutral hydrogen (Henney et al. 2009, eq. A12).
	//CANNOT FIND THIS.
	rate += collisionalExcitationHI(nH, HIIFRAC, T);
	//Collisional ionization equilibrium-cooling curve (Henney et al. 2009, eq. A13).
	if (T > ciec_minT) {
		double cie_rate = ciec*HIIFRAC*HIIFRAC*nH*nH*z0*std::exp(-0.63*std::log(T))*(1.0-std::exp(-std::pow(1.0e-5*T, 1.63)));
		double smoothing = std::min(1.0, (T-5.0e4)/(2.0e4)); //linear smoothing spread over 20000K (PION: cooling.cc).
		rate += cie_rate*smoothing;
	}
	//Neutral and molecular cooling from cloudy models (Henney et al. 2009, eq. A14).
	double T0 = 70.0 + 220.0*std::pow(nH/n0, 0.2);
	rate += nmc*std::pow(nH, 1.6)*std::pow(T, 0.5)*std::exp(-T0/T);
	/*
	std::cout << "imlc = " << imlc*z0*ne*ne*std::exp(-T1/T - (T2/T)*(T2/T)) << "\n";
	std::cout << "nmlc = " << nmlc*z0*ne*nn*std::exp(-T3/T - (T4/T)*(T4/T)) << "\n";
	std::cout << "fffb = " << recombinationHII(nH, HIIFRAC, T) << "\n";
	std::cout << "cxhi = " << collisionalExcitationHI(nH, HIIFRAC, T) << "\n";
	std::cout << "ciec = " << T > ciec_minT ? ciec*HIIFRAC*HIIFRAC*nH*nH*z0*std::exp(-0.63*std::log(T))*(1.0-std::exp(-std::pow(1.0e-5*T, 1.63)))*std::min(1.0, (T-5.0e4)/(2.0e4)) : 0 << "\n";
	std::cout << "nmc = " << nmc*std::pow(nH, 1.6)*std::pow(T, 0.5)*std::exp(-T0/T) << "\n";
	*/
	return rate;
}

/**
 * @brief Calculates the heating rate of gas in a grid cell due to atomic processes.
 * Heating due to ionizing EUV photons; absorption of FUV radiation by dust grains;
 * hard X-rays deep inside the PDR; stellar radiation reprocessed by dense gas (>10^4cm-3)
 * and absorbed by dust; and cosmic ray particles.
 * @param nH Hydrogen number density.
 * @param T Gas temperature.
 * @param F_FUV Flux of far ultra-violet radiation.
 * @param tau Total column density of hydrogen.
 * @return Heating rate [ergs.s-1.cm-3].
 */
double Thermodynamics::heatingProcesses(const double& nH, const double& HIIFRAC, const double& T, const double& Api, const double& F_FUV, const double& tau) {
	double Av_FUV = 1.086*sigmaV*tau; //!< Visual band optical extinction in magnitudes.
	double rate = 0.0;
	//EUV heating (Henney et al. 2009, eq. A2)
	rate += nH*(1.0-HIIFRAC)*Api*excessEnergy;
	//FUV heating (Henney et al. 2009, eq. A3)
	rate += fuvh_a*nH*F_FUV*std::exp(-1.9*Av_FUV)/(fuvh_b +fuvh_c*F_FUV*std::exp(-1.9*Av_FUV)/nH);
	//IR Heating (Henney et al. 2009, eq. A6)
	rate += irh_a*nH*F_FUV*std::exp(-0.05*Av_FUV)*exp(-2.0*log(1.0+irh_b/nH));
	//Cosmic Ray Heating (Henney et al. 2009, eq. A7)
	//Hack: Increasing this by 10X to compensate for no X-ray heating.
	rate += crh*nH;

	return rate;
}

/**
 * @brief Calculates the resultant heating rate of gas in a grid cell due to atomic processes.
 * @param nH Hydrogen number density.
 * @param HIIFRAC Ionized hydrogen fraction.
 * @param T Gas temperature.
 * @param F_FUV Flux of far ultra-violet radiation.
 * @param tau Total column density of hydrogen.
 * @return Resultant heating rate [ergs.s-1.cm-3].
 */
double Thermodynamics::heating(const double& nH, const double& HIIFRAC, const double& T, const double& Api, const double& F_FUV, const double& tau) {
	double rate = 0;
	rate = heatingProcesses(nH, HIIFRAC, T, Api, F_FUV, tau)-coolingProcesses(nH, HIIFRAC, T);
	double T_inf  = 100.0, T_soft = 300.0;
	if (T < T_soft && rate < 0.0) {
		rate = std::min(0.0, rate*(T-T_inf)/(T_soft-T_inf)); //"Soft landing" to equilibrium neutral gas temperature
	}
	return rate;
}

void Thermodynamics::calcHeatFlux(const Radiation& rad) {
	int starID = 0;
	for (GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next) {
		double nH = cptr->Q[iden] / rad.rparams.H_MASS;
		double HIIFRAC = cptr->Q[ihii];
		double T = rad.temperature(*cptr);
		double Api = 0;
		double rsqrd = 0;
		double F_FUV = 0;
		if (rad.stars.size() != 0) {
			Api = rad.photoionisationRate(HIIFRAC, cptr->R[itau], cptr->R[idtau], *cptr, rad.stars[starID]);
			for (int id = 0; id < grid.gparams->ND; ++id)
				rsqrd += cptr->xc[id]-rad.stars[starID].x[id]*grid.dx[id]*grid.dx[id];
			F_FUV = fluxFUV(0.5*rad.stars[0].photonRate, rsqrd);
		}
		double tau = cptr->R[icolden];
		cptr->R[iheat] = heating(nH, HIIFRAC, T, Api, F_FUV, tau);
	}
}

double Thermodynamics::calcTimeStep(const double& dt_dyn) {
	double dt = dt_dyn;
	for (GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next) {
		if (cptr->R[iheat] != 0) {
			dt = std::min(dt, std::abs(0.01*cptr->U[ipre]/cptr->R[iheat]));
		}
	}
	return dt;
}

double Thermodynamics::updateSrcTerms() {
	double udot = 0;
	for (GridCell* cptr = grid.fcell; cptr != NULL; cptr = cptr->next) {
		cptr->UDOT[ipre] += cptr->R[iheat];
		udot = std::abs(cptr->R[iheat]/cptr->U[ipre]) > std::abs(udot) ? cptr->R[iheat]/cptr->U[ipre] : udot;
	}
	return udot;
}
