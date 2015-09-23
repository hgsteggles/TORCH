#include "Radiation.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>
#include <iostream>

#include "Fluid/Fluid.hpp"
#include "Fluid/Grid.hpp"
#include "Fluid/GridCell.hpp"
#include "Fluid/Star.hpp"
#include "IO/Logger.hpp"
#include "MPI/MPI_Wrapper.hpp"
#include "Torch/Converter.hpp"
#include "Torch/Parameters.hpp"

Radiation::Radiation()
: Integrator("Radiation")
{ }

void Radiation::initialise(std::shared_ptr<Constants> c, RadiationParameters rp) {
	m_consts = std::move(c);

	K1 = rp.K1;
	K2 = rp.K2;
	K3 = rp.K3;
	K4 = rp.K4;
	THI = rp.THI;
	THII = rp.THII;
	m_alphaB = rp.alphaB;
	collisions_on = rp.collisions_on;
	heatingAmplification = rp.heatingAmplification;
	massFractionH = rp.massFractionH;
	minX = rp.minX;
	photoIonCrossSection = rp.photoIonCrossSection;
	tau0 = rp.tau0;

	if (rp.coupling.compare("neq") == 0)
		coupling = Coupling::NON_EQUILIBRIUM;
	else if (rp.coupling.compare("tti") == 0)
		coupling = Coupling::TWO_TEMP_ISOTHERMAL;
	else
		coupling = Coupling::OFF;

	if (rp.scheme.compare("implicit2") == 0)
		scheme = Scheme::IMPLICIT2;
	else if (rp.scheme.compare("explicit") == 0)
		scheme = Scheme::EXPLICIT;
	else
		scheme = Scheme::IMPLICIT;

	initRecombinationHummer(m_consts->converter);
}

void Radiation::integrate(double dt, Fluid& fluid) const {
	if (scheme == Scheme::IMPLICIT)
		transferRadiation(dt, fluid);
	else
		transferRadiation2(dt, fluid);
}

int Radiation::getRayPlane(Vec3& xc, Vec3& xs) const {
	int plane = 0;
	double dxcs_max = 0;
	for (int i = 0; i < m_consts->nd; ++i) {
		double dxcs = fabs(xc[i] - xs[i]);
		if (dxcs > dxcs_max) {
			dxcs_max = dxcs;
			plane = i;
		}
	}
	return plane;
}

void Radiation::initField(Fluid& fluid) const {
	static bool time = 0;
	Grid& grid = fluid.getGrid();

	grid.calculateNearestNeighbours(fluid.getStar().xc);
	if (time == 0) {
		for (GridCell& cell : grid.getIterable("GridCells")){
			cell.ds = cellPathLength(cell.xc, fluid.getStar().xc, grid.dx);
			double r_sqrd = 0;
			for(int i = 0; i < m_consts->nd; ++i)
				r_sqrd += (cell.xc[i] - fluid.getStar().xc[i])*(cell.xc[i] - fluid.getStar().xc[i])*grid.dx[i]*grid.dx[i];
			cell.shellVol = shellVolume(cell.ds, r_sqrd);
			double nH = massFractionH*cell.Q[UID::DEN]/m_consts->hydrogenMass;
			cell.R[RID::DTAU] = calc_dtau((1.0 - cell.Q[UID::HII])*nH, cell.ds);
			cell.R[RID::DTAU_A] = calc_dtau((1.0 - cell.R[RID::HII_A])*nH, cell.ds);
			++time;
		}
	}
}

double Radiation::calc_dtau(double nHI, double ds) const {
	return nHI*photoIonCrossSection*ds;
}

double Radiation::cellPathLength(Vec3& xc, Vec3& sc, Vec3& dx) const {
	double denom;
	double d[3] = {0.0, 0.0, 0.0};
	for(int i = 0; i < m_consts->nd; ++i)
		d[i] = std::abs(xc[i] - sc[i]);
	if(d[2] >= d[1] && d[2] >= d[0])
		denom = d[2];
	else if(d[1] > d[2] && d[1] >= d[0])
		denom = d[1];
	else
		denom = d[0];
	if (denom != 0)
		return std::sqrt(d[0]*d[0]*dx[0]*dx[0]+d[1]*d[1]*dx[1]*dx[1]+d[2]*d[2]*dx[2]*dx[2])/denom;
	else
		return 0.0;
}

double Radiation::shellVolume(double ds, double r_sqrd) const {
	//r_min = sqrt(r_sqrd) - 0.5*ds;
	//r_max = sqrt(r_sqrd) + 0.5*ds;
	//if(ND == 1)
	//	return dx*dx*dx;
	//else if(ND == 2)
	//	if(r_sqrd <= 0.1)
	//		return PI*(r_max*r_max - r_min*r_min)*dx;
	//	else
	//		return 2.0*PI*sqrt(r_sqrd)*ds*dx;
	//else if(ND == 3)
	if(r_sqrd <= 100.0*ds*ds) {
		double r = std::sqrt(r_sqrd);
		double r_min = r - 0.5*ds;
		double r_max = r + 0.5*ds;
		return 4.0*m_consts->pi*(r_max*r_max*r_max - r_min*r_min*r_min)/3.0;
	}
	else
		return 4.0*m_consts->pi*r_sqrd*ds;
}

void Radiation::doric(double dt, double& HII_avg, double& HII, double Api, double nHII_aB, double nHII_Aci) const {
	double xeq, inv_ti;
	double epsilon = 1.0e-8;
	//double frac_avg_old;
	double HII_old = HII;
	//frac_avg_old = frac_avg;
	if (nHII_aB == 0.0)
		xeq = 1.0;
	else
		xeq = (Api + nHII_Aci)/(Api + nHII_Aci + nHII_aB);
	inv_ti = Api + nHII_Aci + nHII_aB;
	//if (dt*inv_ti == 0) {
	if (dt*inv_ti < 1.0e-8) {
		/*
		if (inv_ti > 1.0e-7) {
			double HII_exact = xeq + (HII_old-xeq)*exp(-dt*inv_ti);
			double HII_avg_exact = HII_avg = xeq + (HII_old-xeq)*(1.0-exp(-dt*inv_ti))/(dt*inv_ti);
			double HII_error = std::abs(HII_old - HII_exact);
			double HII_avg_error = std::abs(HII_old - HII_avg_exact);
			std::cerr << HII_error << " " << HII_avg_error << '\n';
		}
		 */
		HII = HII_old;
		HII = std::max(std::min(HII, 1.0), 0.0);
		if (1.0 - HII < epsilon)
			HII = 1.0 - epsilon;
		HII_avg = HII;
	}
	else {
		double exp_mdt_inv_ti = std::exp(-dt*inv_ti);
		HII = xeq + (HII_old-xeq)*exp_mdt_inv_ti;
		HII = std::max(std::min(HII, 1.0), 0.0);
		if (1.0 - HII < epsilon)
			HII = 1.0 - epsilon;
		HII_avg = xeq + (HII_old-xeq)*(1.0-exp_mdt_inv_ti)/(dt*inv_ti);
		HII_avg = std::max(std::min(HII_avg, 1.0), 0.0);
		if (1.0 - HII_avg < epsilon)
			HII_avg = 1.0 - epsilon;
	}
}

/**
 * @brief Cubic spline fit for Hummer (1994) HII recombination cooling rate data.
 */
void Radiation::initRecombinationHummer(const Converter& converter) {
	double alphab[31] = { 9.283e-11, 8.823e-11, 8.361e-11, 7.898e-11, 7.435e-11,
			6.973e-11, 6.512e-11, 6.054e-11, 5.599e-11, 5.147e-11, 4.700e-11,
			4.258e-11, 3.823e-11, 3.397e-11, 2.983e-11, 2.584e-11, 2.204e-11,
			1.847e-11, 1.520e-11, 1.226e-11, 9.696e-12, 7.514e-12, 5.710e-12,
			4.257e-12, 3.117e-12, 2.244e-12, 1.590e-12, 1.110e-12, 7.642e-13,
			5.199e-13, 3.498e-13 };
	double coolb[31]  = { 8.287e-11, 7.821e-11, 7.356e-11, 6.892e-11, 6.430e-11,
			5.971e-11, 5.515e-11, 5.062e-11, 4.614e-11, 4.170e-11, 3.734e-11,
			3.306e-11, 2.888e-11, 2.484e-11, 2.098e-11, 1.736e-11, 1.402e-11,
			1.103e-11, 8.442e-12, 6.279e-12, 4.539e-12, 3.192e-12, 2.185e-12,
			1.458e-12, 9.484e-13, 6.023e-13, 3.738e-13, 2.268e-13, 1.348e-13,
			7.859e-14, 4.499e-14 };

	std::vector<std::pair<double, double>> recomb, cool;
	for (int i = 0; i < 31; i++) {
		double T = std::exp(std::log(10.0)*(1.0 + 0.2*static_cast<double>(i)));
		double sqrt_T = std::sqrt(T);
		recomb.push_back(std::make_pair(T, converter.toCodeUnits(alphab[i]/sqrt_T, 0, 3, -1)));
		cool.push_back(std::make_pair(T, converter.toCodeUnits(coolb[i]/sqrt_T, 0, 3, -1)));
	}
	m_recombinationHII_RecombRates = std::unique_ptr<LinearSplineData>(new LinearSplineData(recomb));
	m_recombinationHII_CoolingRates = std::unique_ptr<LinearSplineData>(new LinearSplineData(cool));
}

/**
 * Cubic spline interpolation of the recombination rate coefficient of ionised hydrogen.
 * @param T Gas temperature.
 * @return Recombination cooling rate of HII.
 */
double Radiation::recombinationRateCoefficient(double T) const {
	//return alphaB*pow(T/10000, -0.7); // cm3.s-1
	//return alphaB; // cm3.s-1
	return m_recombinationHII_RecombRates->interpolate(T);
}

/**
 * Cubic spline interpolation of the recombination cooling rate of HII.
 * @param T Gas temperature.
 * @return Recombination cooling rate of HII.
 */
double Radiation::recombinationCoolingRate(double nH, double HIIFRAC, double T) const {
	if (T == 0)
		std::cout << "BOB" << std::endl;
	double rate = m_recombinationHII_CoolingRates->interpolate(T);
	return HIIFRAC*HIIFRAC*nH*nH*m_consts->boltzmannConst*T*rate;
}

/**
 * @brief Calculates the collisional ionisation rate of hydrogen.
 * Uses the fitted formula from Voronov (1997,ADANDT,65,1).
 * @param T Gas temperature.
 * @return Recombination cooling rate of HII.
 */
double Radiation::collisionalIonisationRate(double T) const {
	if (!collisions_on || T < 5.0e3)
		return 0.0;
	double U = m_consts->voronov_DE/T;
	return m_consts->voronov_A*(1.+m_consts->voronov_P*std::sqrt(U))*std::exp(m_consts->voronov_K*log(U) - U)/(m_consts->voronov_X + U); //cm^3/s
}

double Radiation::photoionisationRate(double nHI, double T, double delT, double shellVol, double photonRate) const {
	if (nHI == 0 || shellVol == 0)
		return 0.0;
	else {
		double emT = T > 200 ? 0 : std::exp(-T);
		double emdT = delT > 200 ? 0 : std::exp(-delT);
		if (delT < 1.0e-6)
			return photonRate*emT*delT/(nHI*shellVol);
		return photonRate*emT*(1.0-emdT)/(nHI*shellVol);
	}
}

double Radiation::HIIfracRate(double A_pi, double A_ci, double A_rr, double nH, double frac) const {
	return (1.0-frac)*(A_pi + frac*nH*A_ci) - frac*frac*nH*A_rr;
}

void Radiation::preTimeStepCalculations(Fluid& fluid) const {
	Grid& grid = fluid.getGrid();
	for (int cellID : grid.getOrderedIndices("CausalNonWind")) {
		GridCell& cell = grid.getCell(cellID);

		double n_H = (massFractionH*cell.Q[UID::DEN]/m_consts->hydrogenMass);
		double excessEnergy = fluid.getStar().photonEnergy - m_consts->rydbergEnergy;
		double T = fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]);
		double nHI = (1.0-cell.Q[UID::HII])*n_H;
		double A_pi = photoionisationRate(nHI, cell.R[RID::TAU], cell.R[RID::DTAU], cell.shellVol, fluid.getStar().photonRate);
		double photoion = n_H*(1.0-cell.Q[UID::HII])*A_pi*excessEnergy;
		double recombination = recombinationCoolingRate(n_H, cell.Q[UID::HII], T);
		double collisions = cell.Q[UID::HII]*(1.0-cell.Q[UID::HII])*n_H*n_H*collisionalIonisationRate(T);
		cell.H[HID::RHII] = -recombination;
		double rate = photoion - recombination - collisions;
		if (T < cell.T_min + 200 && rate < 0.0)
			rate = std::min(0.0, rate*(T-cell.T_min)/200); //"Soft landing" to equilibrium neutral gas temperature
		rate *= heatingAmplification;
		cell.R[RID::HEAT] = rate;

		if (rate != rate) {
			std::stringstream out;
			out << "heating rate in radiation module is NaN\n";
			out << "HII = " << cell.Q[UID::HII] << '\n';
			out << "PRE = " << cell.Q[UID::PRE] << '\n';
			out << "DEN = " << cell.Q[UID::DEN] << '\n';
			out << "T = " << T << '\n';
			out << "n_H = " << n_H << '\n';
			out << "A_pi = " << A_pi << '\n';
			out << "phions = " << photoion << '\n';
			out << "recomb = " << recombination << '\n';
			out << "collis = " << collisions << '\n';
			throw std::runtime_error(out.str());
		}
	}
}

double Radiation::calculateTimeStep(double dt_max, Fluid& fluid) const {
	Grid& grid = fluid.getGrid();
	static bool once = false;
	double dt = dt_max, dtc, dt1, dt2, dt3, dt4;
	if (fluid.getStar().on) {
		for (int cellID : grid.getOrderedIndices("CausalNonWind")) {
			GridCell& cell = grid.getCell(cellID);

			dtc = dt1 = dt2 = dt3 = dt4 = dt_max;
			if (coupling == Coupling::NON_EQUILIBRIUM) {
				if (cell.R[RID::HEAT] != 0)
					dtc = std::abs(cell.U[UID::PRE]/cell.R[RID::HEAT]);
			}
			if (K1 != 0.0) {
				double T = fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]);
				if (!once) {
					dt1 = K1*m_consts->hydrogenMass/(massFractionH*cell.Q[UID::DEN]*recombinationRateCoefficient(T));
					once = true;
				}
				else if (cell.Q[UID::HII] != 0)
					dt1 = K1*m_consts->hydrogenMass/(massFractionH*cell.Q[UID::DEN]*cell.Q[UID::HII]*recombinationRateCoefficient(T));
			}
			if(K2 != 0.0)
				dt2 = dt_max;
			if(K3 != 0.0) {
				double nH = massFractionH*cell.Q[UID::DEN]/m_consts->hydrogenMass;
				double nHI = (1.0-cell.Q[UID::HII])*nH;
				double A_pi = photoionisationRate(nHI, cell.R[RID::TAU], cell.R[RID::DTAU], cell.shellVol, fluid.getStar().photonRate);
				double A_ci = collisionalIonisationRate(fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]));
				double A_rr = recombinationRateCoefficient(fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]));
				double fracRate = HIIfracRate(A_pi, A_ci, A_rr, nH, cell.Q[UID::HII]);
				if (fracRate != 0.0)
					dt3 = K3*std::max(0.05, 1.0 - cell.Q[UID::HII])/std::fabs(fracRate);
			}
			if(K4 != 0.0) {
				double nH = massFractionH*cell.Q[UID::DEN]/m_consts->hydrogenMass;
				double nHI = (1.0-cell.Q[UID::HII])*nH;
				double A_pi = photoionisationRate(nHI, cell.R[RID::TAU], cell.R[RID::DTAU], cell.shellVol, fluid.getStar().photonRate);
				double A_ci = collisionalIonisationRate(fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]));
				double A_rr = recombinationRateCoefficient(fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]));
				double fracRate = HIIfracRate(A_pi, A_ci, A_rr, nH, cell.Q[UID::HII]);
				if (fracRate != 0.0)
					dt4 = K4*1.0/fabs(fracRate); // timestep criterion [Mackey 2012].
			}

			dt = std::min(dt, std::min(dt_max, std::min(dtc, std::min(dt1, std::min(dt2, std::min(dt3, dt4))))));

			if (dt == 0.0 || dt != dt || std::isinf(dt))
				throw std::runtime_error("Radiation::calculateTimeStep(): invalid timestep: " + std::to_string(dt));
		}
	}
	return dt;
}

void Radiation::updateTauSC(bool average, GridCell& cell, Fluid& fluid, double dist2) const {
	Grid& grid = fluid.getGrid();
	if(dist2 > 0.95){
		double tau[4] = {0.0, 0.0, 0.0, 0.0};
		double w_raga[4];
		for(int i = 0; i < 4; ++i) {
			if (grid.cellExists(cell.neighbourIDs[i]))
				tau[i] = grid.getCell(cell.neighbourIDs[i]).R[average ? RID::TAU_A : RID::TAU]+grid.getCell(cell.neighbourIDs[i]).R[average ? RID::DTAU_A : RID::DTAU];
			w_raga[i] = cell.neighbourWeights[i]/std::max(tau0, tau[i]);
		}
		double sum_w = w_raga[0]+w_raga[1]+w_raga[2]+w_raga[3];
		double newtau = 0.0;
		for(int i = 0; i < 4; ++i){
			w_raga[i] = w_raga[i]/sum_w;
			newtau += w_raga[i]*tau[i];
		}
		cell.R[average ? RID::TAU_A : RID::TAU] = newtau;
	}
	else
		cell.R[average ? RID::TAU_A : RID::TAU] = 0;
}

void Radiation::update_HIIfrac(double dt, GridCell& cell, Fluid& fluid) const {
	Grid& grid = fluid.getGrid();
	Star& star = fluid.getStar();

	if(!isStar(cell, fluid.getStar())){
		double n_H = massFractionH*cell.Q[UID::DEN] / m_consts->hydrogenMass;
		double A_pi = 0;
		double HII = cell.Q[UID::HII];
		double HII_avg = HII;
		double T = fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]);
		double alphaB = recombinationRateCoefficient(T);
		double A_ci = collisionalIonisationRate(T);
		if (scheme == Scheme::IMPLICIT || scheme == Scheme::IMPLICIT2) {
			double convergence2 = 1.0e-4;
			double convergence_frac = 1.0e-3;
			double HII_avg_old;
			int niter = 0;
			bool converged = false;
			//HII_avg = cell.R[ihiita];
			double tau_avg = cell.R[RID::TAU_A];
			if (scheme == Scheme::IMPLICIT2) tau_avg = cell.R[RID::TAU];
			while (!converged){
				niter++;
				HII_avg_old = HII_avg;
				HII = cell.Q[UID::HII];
				double dtau_avg = (1.0-HII_avg)*n_H*photoIonCrossSection*cell.ds;
				//double T = temperature(cell.Q[ipre], cell.Q[iden], HII_avg);
				double nHI = (1.0-HII_avg)*n_H;
				A_pi = photoionisationRate(nHI, tau_avg, dtau_avg, cell.shellVol, fluid.getStar().photonRate);
				double nHII_aB = HII_avg*n_H*alphaB;
				double nHII_Aci = HII_avg*n_H*A_ci;

				doric(dt, HII, HII_avg, A_pi, nHII_aB, nHII_Aci);

				if ((fabs((HII_avg-HII_avg_old)/HII_avg) < convergence2 || (HII_avg < convergence_frac) || A_pi == 0))
					converged = true;
				if (niter > 50000) {
					Logger<FileLogPolicy>::Instance().print<SeverityType::WARNING>(
							"Radiation::calculate_HIIfrac: slow iteration: HII_avg = ", HII_avg, ", HII_avg_old = ", HII_avg_old, ", HII = ", HII
					);
				}
				if (niter > 50020 || HII != HII) {
					std::stringstream out;
					out << "Radiation::calculate_HIIfrac: implicit method not converging.\n";
					out << "temperature = " << fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]) << '\n';
					out << "A_pi = " << A_pi << '\n';
					out << "alphaB = " << alphaB << '\n';
					out << "hii_dot = " << HIIfracRate(A_pi, A_ci, alphaB, n_H, cell.Q[UID::HII]) << '\n';
					out << printInfo();
					out << cell.printInfo();
					if (cell.R[RID::TAU] != cell.R[RID::TAU]) {
						out << "Radiation::calculate_HIIfrac(): tau is NaN.\n";
						for (int neighbourID : cell.neighbourIDs) {
							if (grid.cellExists(neighbourID))
								out << grid.getCell(neighbourID).printInfo();
						}
					}
					throw std::runtime_error(out.str());
				}
			}
			cell.R[RID::HII_A] = HII_avg;
		}
		else if (scheme == Scheme::EXPLICIT){
			double HII_dummy, tau, dtau;
			tau = cell.R[RID::TAU];
			dtau = cell.R[RID::DTAU];
			double nHI = (1.0-HII)*(massFractionH*cell.Q[UID::DEN]/m_consts->hydrogenMass);
			A_pi = photoionisationRate(nHI, tau, dtau, cell.shellVol, fluid.getStar().photonRate);
			double nHII_aB = HII_avg*n_H*alphaB;
			double nHII_Aci = HII_avg*n_H*A_ci;
			//set_HIIfrac(HII+dt*HIIfracDot(A_pi, HII) );
			HII_dummy = HII;
			doric(dt, HII, HII_dummy, A_pi, nHII_aB, nHII_Aci);
		}

		/** Photoionisation heating and cooling from collisional ionisations and recombinations **/
		/*
		double excessEnergy = fluid.getStar().sparams->photonEnergy - consts->rydbergEnergy;
		//double T = temperature(cell.Q[ipre], cell.Q[iden], HII_avg);
		double pion_heating = n_H*(1.0-HII_avg)*A_pi*excessEnergy;
		double recomb_cooling = recombinationCoolingRate(n_H, HII_avg, T);
		double coll_cooling = HII_avg*(1.0-HII_avg)*n_H*n_H*A_ci;

		cell.R[iheat] = pion_heating - recomb_cooling - coll_cooling;
		if (T < 300 && cell.R[iheat] < 0.0)
			cell.R[iheat] = std::min(0.0, cell.R[iheat]*(T-100)/(300-100)); //"Soft landing" to equilibrium neutral gas temperature
		cell.R[iheat] *= heatingAmplification;
		*/
		cell.Q[UID::HII] = HII;
		if (HII > 0.001)
			cell.Q[UID::ADV] = 1;

		/*
		cell.H[0] = A_pi;
		cell.H[1] = alphaB;
		cell.H[2] = cell.R[RID::DTAU];
		cell.H[3] = cell.R[RID::TAU];
		cell.H[4] = cell.R[RID::DTAU_A];
		cell.H[5] = cell.R[RID::TAU_A];
		*/
	}
	else {
		cell.R[RID::HEAT] = 0;
		cell.Q[UID::HII] = 1;
		cell.R[RID::HII_A] = 1;
	}
}

void Radiation::rayTrace(Fluid& fluid) const {
	Grid& grid = fluid.getGrid();
	if (fluid.getStar().on) {
		for (int cellID : grid.getOrderedIndices("CausalWind")) {
			GridCell& cell = grid.getCell(cellID);

			cell.R[RID::TAU] = 0;
			cell.R[RID::TAU_A] = 0;
			cell.R[RID::DTAU] = 0;
			cell.R[RID::DTAU_A] = 0;
			cell.Q[UID::HII] = 1;
			cell.R[RID::HII_A] = 1;
		}
		bool average = true;
		/** Causally loop over cells in grid */
		for (int cellID : grid.getOrderedIndices("CausalNonWind")) {
			GridCell& cell = grid.getCell(cellID);

			double dist2 = 0;
			for (int i = 0; i < m_consts->nd; ++i)
				dist2 += (cell.xc[i] - fluid.getStar().xc[i])*(cell.xc[i] - fluid.getStar().xc[i]);
			/** Calculate column densities */
			updateTauSC(average==false, cell, fluid, dist2);
			double nH = massFractionH*cell.Q[UID::DEN]/m_consts->hydrogenMass;
			cell.R[RID::DTAU] = calc_dtau((1.0 - cell.Q[UID::HII])*nH, cell.ds);
			cell.R[RID::DTAU_A] = calc_dtau((1.0 - cell.R[RID::HII_A])*nH, cell.ds);
		}
	}
}

void Radiation::transferRadiation2(double dt, Fluid& fluid) const {
	MPIW& mpihandler = MPIW::Instance();
	Grid& grid = fluid.getGrid();
	Star& star = fluid.getStar();

	if (star.on) {
		PartitionManager& partition = grid.getPartitionManager();
		/** Wait for column densities */
		if (star.core != Star::Location::HERE) {
			int source = star.core == Star::Location::LEFT ? mpihandler.rank - 1 : mpihandler.rank + 1;
			partition.recvData(source, SendID::RADIATION_MSG);

			for (GridCell& ghost : grid.getIterable(star.core == Star::Location::LEFT ? "LeftPartitionCells" : "RightPartitionCells")) {
				ghost.R[RID::DTAU] = partition.getRecvItem();
				ghost.R[RID::TAU] = partition.getRecvItem();
			}
		}
		/** Calculate column densities */
		rayTrace(fluid);
		/** Send column densities to processor on left */
		if (!(mpihandler.getRank() == 0 || star.core == Star::Location::LEFT)) {
			for (GridCell& ghost : grid.getIterable("LeftPartitionCells")) {
				GridCell& cell = grid.right(0, ghost);
				partition.addSendItem(cell.R[RID::DTAU]);
				partition.addSendItem(cell.R[RID::TAU]);
			}
			int destination = mpihandler.rank - 1;
			partition.sendData(destination, SendID::RADIATION_MSG);
		}
		/** Send column densities to processor on right */
		if (!(mpihandler.getRank() == mpihandler.nProcessors()-1 || star.core == Star::Location::RIGHT)) {
			for (GridCell& ghost : grid.getIterable("RightPartitionCells")) {
				GridCell& cell = grid.left(0, ghost);
				partition.addSendItem(cell.R[RID::DTAU]);
				partition.addSendItem(cell.R[RID::TAU]);
			}
			int destination = mpihandler.rank + 1;
			partition.sendData(destination, SendID::RADIATION_MSG);
		}
		for (int cellID : grid.getOrderedIndices("CausalNonWind")) {
			GridCell& cell = grid.getCell(cellID);
			update_HIIfrac(dt, cell, fluid);
		}
	}
	else {
		double HII_dummy = 0, HII;
		for (GridCell& cell : grid.getIterable("GridCells")) {
			HII = cell.Q[UID::HII];
			HII_dummy = HII;
			double n_H = massFractionH*cell.Q[UID::DEN] / m_consts->hydrogenMass;
			double T = fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]);
			double nHII_aB = HII*n_H*recombinationRateCoefficient(T);
			double nHII_Aci = HII*n_H*collisionalIonisationRate(T);
			doric(dt, HII, HII_dummy, 0, nHII_aB, nHII_Aci);
			cell.Q[UID::HII] = HII;
		}
	}
}

void Radiation::transferRadiation(double dt, Fluid& fluid) const {
	MPIW& mpihandler = MPIW::Instance();
	Grid& grid = fluid.getGrid();
	Star& star = fluid.getStar();

	if (fluid.getStar().on) {
		PartitionManager& partition = grid.getPartitionManager();
		partition.resetBuffer();
		/** Wait for column densities */
		if (fluid.getStar().core != Star::Location::HERE) {
			int source = star.core == Star::Location::LEFT ? mpihandler.rank - 1 : mpihandler.rank + 1;
			partition.recvData(source, SendID::RADIATION_MSG);

			for (GridCell& ghost : grid.getIterable(star.core == Star::Location::LEFT ? "LeftPartitionCells" : "RightPartitionCells")) {
				ghost.R[RID::DTAU] = partition.getRecvItem();
				ghost.R[RID::TAU] = partition.getRecvItem();
				ghost.R[RID::DTAU_A] = partition.getRecvItem();
				ghost.R[RID::TAU_A] = partition.getRecvItem();
			}
		}
		/** Causal ray tracing and integrating for HII fraction */
		for (int cellID : grid.getOrderedIndices("CausalWind")) {
			GridCell& cell = grid.getCell(cellID);
			cell.R[RID::TAU] = 0;
			cell.R[RID::TAU_A] = 0;
			cell.R[RID::DTAU] = 0;
			cell.R[RID::DTAU_A] = 0;
			cell.Q[UID::HII] = 1;
			cell.R[RID::HII_A] = 1;
		}
		bool average = true;
		for (int cellID : grid.getOrderedIndices("CausalNonWind")) {
			GridCell& cell = grid.getCell(cellID);

			double dist2 = 0;
			for (int i = 0; i < m_consts->nd; ++i)
				dist2 += (cell.xc[i] - fluid.getStar().xc[i])*(cell.xc[i] - fluid.getStar().xc[i]);
			updateTauSC(average==false, cell, fluid, dist2);
			updateTauSC(average==true, cell, fluid, dist2);
			update_HIIfrac(dt, cell, fluid);
			double nH = massFractionH*cell.Q[UID::DEN]/m_consts->hydrogenMass;
			cell.R[RID::DTAU] = calc_dtau((1.0 - cell.Q[UID::HII])*nH, cell.ds);
			cell.R[RID::DTAU_A] = calc_dtau((1.0 - cell.R[RID::HII_A])*nH, cell.ds);
		}
		/** Send column densities to processor on left */
		if (!(mpihandler.getRank() == 0 || fluid.getStar().core == Star::Location::LEFT)) {
			for (GridCell& ghost : grid.getIterable("LeftPartitionCells")) {
				GridCell& cell = grid.right(0, ghost);
				partition.addSendItem(cell.R[RID::DTAU]);
				partition.addSendItem(cell.R[RID::TAU]);
				partition.addSendItem(cell.R[RID::DTAU_A]);
				partition.addSendItem(cell.R[RID::TAU_A]);
			}
			int destination = mpihandler.rank - 1;
			partition.sendData(destination, SendID::RADIATION_MSG);
		}
		/** Send column densities to processor on right */
		if (!(mpihandler.getRank() == mpihandler.nProcessors()-1 || fluid.getStar().core == Star::Location::RIGHT)) {
			for (GridCell& ghost : grid.getIterable("RightPartitionCells")) {
				GridCell& cell = grid.left(0, ghost);
				partition.addSendItem(cell.R[RID::DTAU]);
				partition.addSendItem(cell.R[RID::TAU]);
				partition.addSendItem(cell.R[RID::DTAU_A]);
				partition.addSendItem(cell.R[RID::TAU_A]);
			}
			int destination = mpihandler.rank + 1;
			partition.sendData(destination, SendID::RADIATION_MSG);
		}
	}
	else {
		double HII_dummy = 0, HII;
		for (GridCell& cell : grid.getIterable("AllCells")) {
			HII = cell.Q[UID::HII];
			HII_dummy = HII;
			double n_H = massFractionH*cell.Q[UID::DEN] / m_consts->hydrogenMass;
			double T = fluid.calcTemperature(cell.Q[UID::HII], cell.Q[UID::PRE], cell.Q[UID::DEN]);
			double nHII_aB = HII*n_H*recombinationRateCoefficient(T);
			double nHII_Aci = HII*n_H*collisionalIonisationRate(T);
			doric(dt, HII, HII_dummy, 0, nHII_aB, nHII_Aci);
			cell.Q[UID::HII] = HII;
		}
	}
}

void Radiation::updateSourceTerms(double dt, Fluid& fluid) const {
	Grid& grid = fluid.getGrid();
	if (fluid.getStar().on) {
		for (GridCell& cell : grid.getIterable("GridCells")) {
			if (coupling == Coupling::TWO_TEMP_ISOTHERMAL) {
				double HII_new = cell.Q[UID::HII];
				double mu_inv_new = massFractionH*(HII_new + 1.0) + (1.0 - massFractionH)*0.25;
				double T_new = THI + (THII-THI)*HII_new;
				double pre_new = m_consts->specificGasConstant*mu_inv_new*cell.Q[UID::DEN]*T_new;

				double ke = 0.0;
				for (int dim = 0; dim < m_consts->nd; ++dim)
					ke += 0.5*cell.U[UID::VEL+dim]*cell.U[UID::VEL+dim]/cell.U[UID::DEN];


				double E_new = pre_new/(cell.heatCapacityRatio-1.0) + ke;

				cell.UDOT[UID::PRE] += (E_new-cell.U[UID::PRE])/dt;
			}
			else if (coupling == Coupling::NON_EQUILIBRIUM) {
				cell.UDOT[UID::PRE] += cell.R[RID::HEAT];
			}

			cell.UDOT[UID::HII] += (cell.Q[UID::HII]*cell.Q[UID::DEN] - cell.U[UID::HII])/dt;
			cell.UDOT[UID::ADV] += (cell.Q[UID::ADV]*cell.Q[UID::DEN] - cell.U[UID::ADV])/dt;
		}
	}
}

bool Radiation::isStar(const GridCell& cell, const Star& star) const {
	return ((int)cell.xc[0] == (int)star.xc[0] && (int)cell.xc[1] == (int)star.xc[1] && (int)cell.xc[2] == star.xc[2]);
}

std::string Radiation::printInfo() const {
	std::stringstream out;

	out << "K1 = " << K1 << "\n";
	out << "K2 = " << K2 << "\n";
	out << "K3 = " << K3 << "\n";
	out << "K4 = " << K4 << "\n";
	out << "photoIonCrossSection = " << photoIonCrossSection << "\n";
	out << "alphaB = " << m_alphaB << "\n";
	out << "tau0 = " << tau0 << "\n";
	out << "minX = " << minX << "\n";
	out << "THI = " << THI << "\n";
	out << "THII = " << THII << "\n";
	out << "scheme = " << (int)scheme << "\n";
	out << "massFractionH = " << massFractionH << "\n";
	out << "heatingAmplification = " << heatingAmplification << "\n";
	out << "collisions_on = " << collisions_on << "\n";
	out << "coupling = " << (int)coupling << "\n";

	return out.str();
}
