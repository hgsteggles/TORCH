#include "GridCell.hpp"

#include <cmath>
#include <iostream>
#include <sstream>

/**
 * @brief The default GridJoin constructor.
 * Provides all attributes with safe values.
 */
GridJoin::GridJoin() {
	for(int i = 0; i < UID::N; ++i)
		F[i] = 0;
	s_total++;
}

/**
 * @brief The Gridjoin destructor.
 */
GridJoin::~GridJoin() {
	s_total--;
}

int GridJoin::s_total = 0;

/**
 * @brief Default GridCell constructor.
 * Provides all attributes with safe values.
 */
GridCell::GridCell() {
	for (int dim = 0; dim < 3; ++dim) {
		for (int iu = 0; iu < UID::N; ++iu) {
			QL[dim][iu] = 0;
			QR[dim][iu] = 0;
		}
		GRAV[dim] = 0;
	}
	for (int i = 0; i < UID::N; ++i) {
		UDOT[i] = 0;
		U[i] = 0;
		Q[i] = 0;
		W[i] = 0;
	}
	for (int i = 0; i < RID::N; ++i)
		R[i] = 0;
	for (int i = 0; i < TID::N; ++i)
		T[i] = 0;
	for (int i = 0; i < HID::N; ++i)
		H[i] = 0;
	s_total++;
}

/**
 * @brief GridCell destructor.
 * Does NOT delete objects that this GridCell object points to.
 */
GridCell::~GridCell() {
	s_total--;
}

void GridCell::setSoundSpeed(double a) {
	m_soundSpeed = a;
}

double GridCell::getSoundSpeed() const {
	return m_soundSpeed;
}

std::string GridCell::printCoords() const {
	std::stringstream out;
	out << xc[0] << ", " << xc[1] << ", " << xc[2] << '\n';
	return out.str();
}

std::string GridCell::printInfo() const {
	std::stringstream out;
	out << "xc[0] = " << xc[0] << '\n';
	out << "xc[1] = " << xc[1] << '\n';
	out << "xc[2] = " << xc[2] << '\n';
	out << "den = " << Q[UID::DEN] << '\n';
	out << "pre = " << Q[UID::PRE] << '\n';
	out << "vel+0 = " << Q[UID::VEL+0] << '\n';
	out << "vel+1 = " << Q[UID::VEL+1] << '\n';
	out << "vel+2 = " << Q[UID::VEL+2] << '\n';
	out << "hii = " << Q[UID::HII] << '\n';
	out << "adv = " << Q[UID::ADV] << '\n';

	for (int i = 0; i < 3; ++i)
		out << "l_den[" << i << "] = " << QL[i][UID::DEN] << '\n';
	for (int i = 0; i < 3; ++i)
		out << "l_pre[" << i << "] = " << QL[i][UID::PRE] << '\n';
	for (int j = 0; j < 3; ++j) {
		for (int i = 0; i < 3; ++i)
			out << "l_vel" << j << "[" << i << "] = " << QL[i][UID::VEL+j] << '\n';
	}
	for (int i = 0; i < 3; ++i)
		out << "l_hii[" << i << "] = " << QL[i][UID::HII] << '\n';
	for (int i = 0; i < 3; ++i)
			out << "l_adv[" << i << "] = " << QL[i][UID::ADV] << '\n';

	for (int i = 0; i < 3; ++i)
		out << "r_den[" << i << "] = " << QR[i][UID::DEN] << '\n';
	for (int i = 0; i < 3; ++i)
		out << "r_pre[" << i << "] = " << QR[i][UID::PRE] << '\n';
	for (int j = 0; j < 3; ++j) {
		for (int i = 0; i < 3; ++i)
			out << "r_vel" << j << "[" << i << "] = " << QR[i][UID::VEL+j] << '\n';
	}
	for (int i = 0; i < 3; ++i)
		out << "r_hii[" << i << "] = " << QR[i][UID::HII] << '\n';
	for (int i = 0; i < 3; ++i)
		out << "r_adv[" << i << "] = " << QR[i][UID::ADV] << '\n';

	out << "u_den = " << U[UID::DEN] << '\n';
	out << "u_pre = " << U[UID::PRE] << '\n';
	out << "u_vel+0 = " << U[UID::VEL+0] << '\n';
	out << "u_vel+1 = " << U[UID::VEL+1] << '\n';
	out << "u_vel+2 = " << U[UID::VEL+2] << '\n';
	out << "u_hii = " << U[UID::HII] << '\n';
	out << "u_adv = " << U[UID::ADV] << '\n';
	out << "heatCapacityRatio = " << heatCapacityRatio << '\n';
	out << "tau = " << R[RID::TAU] << '\n';
	out << "tau_a = " << R[RID::TAU_A] << '\n';
	out << "dtau = " << R[RID::DTAU] << '\n';
	out << "dtau_a = " << R[RID::DTAU_A] << '\n';
	out << "NN[0] = " << NN[0] << '\n';
	out << "NN[1] = " << NN[1] << '\n';
	out << "NN[2] = " << NN[2] << '\n';
	out << "NN[3] = " << NN[3] << '\n';
	out << "NN_weights[0] = " << NN_weights[0] << '\n';
	out << "NN_weights[1] = " << NN_weights[1] << '\n';
	out << "NN_weights[2] = " << NN_weights[2] << '\n';
	out << "NN_weights[3] = " << NN_weights[3] << '\n';
	out << "ljoin[0] = " << ljoin[0] << '\n';
	out << "ljoin[1] = " << ljoin[1] << '\n';
	out << "ljoin[2] = " << ljoin[2] << '\n';
	out << "rjoin[0] = " << rjoin[0] << '\n';
	out << "rjoin[1] = " << rjoin[1] << '\n';
	out << "rjoin[2] = " << rjoin[2] << '\n';
	out << "left[0] = " << left[0] << '\n';
	out << "left[1] = " << left[1] << '\n';
	out << "left[2] = " << left[2] << '\n';
	out << "right[0] = " << right[0] << '\n';
	out << "right[1] = " << right[1] << '\n';
	out << "right[2] = " << right[2] << '\n';

	return out.str();
}

int GridCell::s_total = 0;

/**
 * @brief Setter for GridCell::U.
 * @param index
 * @param value
 */
void GridCell::set_U(const int index, const double value) {U[index] = value;}

/**
 * @brief Setter for GridCell:xc.
 * @param x The x grid coordinate.
 * @param y The y grid coordinate.
 * @param z The z grid coordinate.
 */
void GridCell::set_xcs(const double x, const double y, const double z) {
	xc[0] = x;
	xc[1] = y;
	xc[2] = z;
}

/**
 * @brief Getter for GridCell::xc.
 * @param i The grid coordinate to be returned.
 * @return The location of this GridCell object on grid coordinate i.
 */
double GridCell::get_xc(const int index) {return xc[index];}

/**
 * @brief Getter for GridCell::U.
 * @param index The index for the fluid variable to be returned.
 * @return The value of the fluid variable.
 */
double GridCell::get_U(const int index) {return U[index];}

double GridCell::temperature(const double massFracH, const double specGasConst) const {
	double mu_inv = massFracH*(Q[UID::HII] + 1.0) + (1.0 - massFracH)*0.25;
	return (Q[UID::PRE]/Q[UID::DEN])*(1.0/mu_inv)/specGasConst;
}

void UfromQ(FluidArray& u, const FluidArray& q, double gamma, int nd) {
	double ke = 0;
	if (!std::isfinite(q[UID::DEN]) || !std::isfinite(q[UID::PRE]))
		std::cout << "ufromq: " << q[UID::DEN] << '\t' << q[UID::PRE] << std::endl;
	u[UID::DEN] = q[UID::DEN];
	for(int dim = 0; dim < nd; ++dim){
		u[UID::VEL+dim] = q[UID::VEL+dim]*q[UID::DEN];
		ke += 0.5*q[UID::DEN]*q[UID::VEL+dim]*q[UID::VEL+dim];
	}
	u[UID::PRE] = q[UID::PRE]/(gamma - 1.0) + ke;
	u[UID::HII] = q[UID::HII]*q[UID::DEN];
	u[UID::ADV] = q[UID::ADV]*q[UID::DEN];
}

void QfromU(FluidArray& q, const FluidArray& u, double gamma, int nd) {
	q[UID::DEN] = u[UID::DEN];
	double ke = 0;
	for (int dim = 0; dim < nd; ++dim) {
		q[UID::VEL+dim] = u[UID::VEL+dim]/u[UID::DEN];
		ke += 0.5*u[UID::VEL+dim]*u[UID::VEL+dim]/u[UID::DEN];
	}
	q[UID::PRE] = (u[UID::PRE] - ke)*(gamma - 1.0);
	q[UID::HII] = u[UID::HII]/u[UID::DEN];
	q[UID::ADV] = u[UID::ADV]/u[UID::DEN];
}

void FfromU(FluidArray& f, const FluidArray& u, double gamma, int nd, int dim) {
	f[UID::DEN] = u[UID::VEL+dim];
	double ke = 0, pressure;
	for(int id = 0; id < nd; ++id) {
		f[UID::VEL+(dim+id)%nd] = u[UID::VEL+(dim+id)%nd]*u[UID::VEL+dim]/u[UID::DEN];
		ke += 0.5*u[UID::VEL+id]*u[UID::VEL+id]/u[UID::DEN];
	}
	pressure = (u[UID::PRE] - ke)*(gamma - 1.0);
	f[UID::VEL+dim] += pressure;
	f[UID::PRE] = u[UID::VEL+dim]*(u[UID::PRE] + pressure)/u[UID::DEN];
	f[UID::HII] = u[UID::VEL+dim]*u[UID::HII]/u[UID::DEN];
	f[UID::ADV] = u[UID::VEL+dim]*u[UID::ADV]/u[UID::DEN];
}
void FfromQ(FluidArray& f, const FluidArray& q, double gamma, int nd, const int dim) {
	f[UID::DEN] = q[UID::DEN]*q[UID::VEL+dim];
	double ke = 0;
	for(int id = 0; id < nd; ++id) {
		ke += q[UID::VEL+id]*q[UID::VEL+id];
		f[UID::VEL+(dim+id)%nd] = q[UID::DEN]*q[UID::VEL+(dim+id)%nd]*q[UID::VEL+dim];
	}
	ke *= 0.5*q[UID::DEN];
	f[UID::VEL+dim] += q[UID::PRE];
	double g2 = gamma/(gamma - 1.0);
	f[UID::PRE] = q[UID::VEL+dim]*(g2*q[UID::PRE] + ke);
	f[UID::HII] = q[UID::DEN]*q[UID::VEL+dim]*q[UID::HII];
	f[UID::ADV] = q[UID::DEN]*q[UID::VEL+dim]*q[UID::ADV];
}

