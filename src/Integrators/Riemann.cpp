#include "Riemann.hpp"

#include "Eigen/Dense"

#include <cmath>
#include <iostream>

std::string printQ(const FluidArray& Q) {
	std::stringstream out;

	out << "density  = " << Q[UID::DEN] << std::endl;
	out << "pressure = " << Q[UID::PRE] << std::endl;
	out << "hii      = " << Q[UID::HII] << std::endl;
	out << "vel0     = " << Q[UID::VEL+0] << std::endl;
	out << "vel1     = " << Q[UID::VEL+1] << std::endl;
	out << "vel2     = " << Q[UID::VEL+2] << std::endl;

	return out.str();
}

std::pair<double, double> characteristicWaveSpeeds(double a_l2, double a_r2, const FluidArray& Q_l,  const FluidArray& Q_r, const double gamma, const int dim) {
	/**
	 * Einfeldt Estimates
	 * Used in HLLE solver
	 * Toro: Riemann Solvers + Numerical Methods for Fluid Dynamics (pg 328)
	 */
	double sqrtrho_l = std::sqrt(Q_l[UID::DEN]);
	double sqrtrho_r = std::sqrt(Q_r[UID::DEN]);
	double u_tilde = (sqrtrho_l*Q_l[UID::VEL+dim] + sqrtrho_r*Q_r[UID::VEL+dim])/(sqrtrho_l + sqrtrho_r);
	double nu2 = 0.5*sqrtrho_r*sqrtrho_l/((sqrtrho_l + sqrtrho_r)*(sqrtrho_l + sqrtrho_r));
	double dsqrd = (sqrtrho_l*a_l2 + sqrtrho_r*a_r2)/(sqrtrho_l + sqrtrho_r);
	dsqrd += nu2*(Q_r[UID::VEL+dim] - Q_l[UID::VEL+dim])*(Q_r[UID::VEL+dim] - Q_l[UID::VEL+dim]);
	double d = std::sqrt(dsqrd);

	return std::make_pair<double, double>(u_tilde - d, u_tilde + d);
}

RiemannSolver::RiemannSolver(int nd) : m_ND(nd) { }

int RiemannSolver::getNumberDimensions() const {
	return m_ND;
}

HartenLaxLeerContactSolver::HartenLaxLeerContactSolver(int nd) : RiemannSolver(nd) { }

void HartenLaxLeerContactSolver::solve(FluidArray& F, const FluidArray& Q_l, const FluidArray& Q_r, double a_l2, double a_r2, double gamma, int dim) const {
	std::pair<double, double> S = characteristicWaveSpeeds(a_l2, a_r2, Q_l, Q_r, gamma, dim);
	double S_l = S.first;
	double S_r = S.second;
	int nd = getNumberDimensions();
	if (S_l >= 0) {
		FfromQ(F, Q_l, gamma, nd, dim);
	}
	else if (S_r <= 0) {
		FfromQ(F, Q_r, gamma, nd, dim);
	}
	else {
		double S_c = Q_r[UID::PRE]-Q_l[UID::PRE];
		S_c += Q_l[UID::DEN]*Q_l[UID::VEL+dim]*(S_l-Q_l[UID::VEL+dim]);
		S_c -= Q_r[UID::DEN]*Q_r[UID::VEL+dim]*(S_r-Q_r[UID::VEL+dim]);
		S_c /= (Q_l[UID::DEN]*(S_l-Q_l[UID::VEL+dim])-Q_r[UID::DEN]*(S_r-Q_r[UID::VEL+dim]));

		const FluidArray& Q_lr = (S_c >= 0) ? Q_l : Q_r;
		double S_lr = (S_c >= 0) ? S_l : S_r;

		FluidArray U_lr, F_lr;
		for (unsigned int i = 0; i < U_lr.size(); ++i)
			U_lr[i] = F_lr[i] = 0;
		UfromQ(U_lr, Q_lr, gamma, nd);
		FfromQ(F_lr, Q_lr, gamma, nd, dim);

		double A_lr = Q_lr[UID::DEN]*(S_lr-Q_lr[UID::VEL+dim])/(S_lr-S_c);
		FluidArray U_clr;
		U_clr[UID::DEN] = A_lr;
		for (int id = 0; id < nd; ++id)
			U_clr[UID::VEL+id] = A_lr*Q_lr[UID::VEL+id];
		U_clr[UID::VEL+dim] = A_lr*S_c;
		U_clr[UID::PRE] = A_lr*((U_lr[UID::PRE]/Q_lr[UID::DEN]) + (S_c-Q_lr[UID::VEL+dim])*(S_c+Q_lr[UID::PRE]/(Q_lr[UID::DEN]*(S_lr-Q_lr[UID::VEL+dim]))));
		U_clr[UID::HII] = A_lr*Q_lr[UID::HII];

		for (int i = 0; i < UID::N; ++i)
			F[i] = F_lr[i] + S_lr*(U_clr[i] - U_lr[i]);
		if (nd < 3)
			F[UID::VEL+2] = 0;
		if (nd < 2)
			F[UID::VEL+1] = 0;
	}
	if ( F[UID::HII] != F[UID::HII] ) {
		std::stringstream out;
		out << "HartenLaxLeerContactSolver::solve: HII flux is NaN\n";
		for (int iu = 0; iu < UID::N; ++iu)
			out << "F[" << iu << "] = " << F[iu] << '\n';
		out << "Q_l is: \n" << printQ(Q_l) << "Q_r is: \n" << printQ(Q_r);

		throw std::runtime_error(out.str());
	}
}

HartenLaxLeerSolver::HartenLaxLeerSolver(int nd) : RiemannSolver(nd) { }

void HartenLaxLeerSolver::solve(FluidArray& F, const FluidArray& Q_l, const FluidArray& Q_r, double a_l2, double a_r2, double gamma, int dim) const {
	std::pair<double, double> S = characteristicWaveSpeeds(a_l2, a_r2, Q_l, Q_r, gamma, dim);
	double S_l = S.first;
	double S_r = S.second;
	if (S_r <= 0) {
		FfromQ(F, Q_r, gamma, getNumberDimensions(), dim);
	}
	else if (S_l >= 0) {
		FfromQ(F, Q_l, gamma, getNumberDimensions(), dim);
	}
	else {
		FluidArray U_l, U_r, F_l, F_r;
		for(int i = 0; i < UID::N; ++i)
			U_l[i] = U_r[i] = F_l[i] = F_r[i] = 0;
		UfromQ(U_l, Q_l, gamma, getNumberDimensions());
		UfromQ(U_r, Q_r, gamma, getNumberDimensions());
		FfromQ(F_l, Q_l, gamma, getNumberDimensions(), dim);
		FfromQ(F_r, Q_r, gamma, getNumberDimensions(), dim);
		for (int i = 0; i < UID::N; ++i)
			F[i] = (S_r*F_l[i] - S_l*F_r[i] + S_l*S_r*(U_r[i]-U_l[i]))/(S_r-S_l);
	}
	if ( F[UID::HII] != F[UID::HII] ) {
		std::stringstream out;
		out << "HartenLaxLeerSolver::solve: HII flux is NaN.\n";
		for (int iu = 0; iu < UID::N; ++iu)
			out << "F[" << iu << "] = " << F[iu] << '\n';
		out << "Q_l is: \n";
		out << printQ(Q_l);
		out << "Q_r is: \n";
		out << printQ(Q_r);
		throw std::runtime_error(out.str());
	}
}

Eigen::Matrix<double, 3, 3> getRotationMatrix(Eigen::Matrix<double, 3, 1> l, double cos, double sin) {
	Eigen::Matrix<double, 3, 3> R;
	R << l(0)*l(0) + (l(1)*l(1) + l(2)*l(2))*cos, l(0)*l(1)*(1.0-cos) - l(2)*sin         , l(0)*l(2)*(1.0-cos) + l(1)*sin         ,
		 l(0)*l(1)*(1.0-cos) + l(2)*sin         , l(1)*l(1) + (l(0)*l(0) + l(2)*l(2))*cos, l(1)*l(2)*(1.0-cos) - l(0)*sin         ,
		 l(0)*l(2)*(1.0-cos) - l(1)*sin         , l(1)*l(2)*(1.0-cos) + l(0)*sin         , l(2)*l(2) + (l(0)*l(0) + l(1)*l(1))*cos;
	return R;
}

Eigen::Matrix<double, 3, 3> getInverseRotationMatrix(Eigen::Matrix<double, 3, 1> l, double cos, double sin) {
	return getRotationMatrix(l, cos, -sin);
}

void rotate(Eigen::Matrix<double, 3, 3> R, std::array<double, UID::N>& v) {
	Eigen::Matrix<double, 3, 1> vel(v[UID::VEL+0], v[UID::VEL+1], v[UID::VEL+2]);
	vel = R*vel;
	for (int i = 0; i < 3; ++i)
		v[UID::VEL+i] = vel(i);
}

RotatedHartenLaxLeerSolver::RotatedHartenLaxLeerSolver(int nd)
	: RiemannSolver(nd)
	, m_hllc(nd)
	, m_hll(nd)
{

}

void RotatedHartenLaxLeerSolver::solve(FluidArray& F, const FluidArray& Q_l, const FluidArray& Q_r, double a_l2, double a_r2, double gamma, int dim) const {
	bool debugger = false;
	bool pureHLLC = false, pureHLL = false;

	Eigen::Matrix<double, 3, 1> d(dim==0 ? 1 : 0, dim==1 ? 1 : 0, dim==2 ? 1 : 0);
	Eigen::Matrix<double, 3, 1> n1(Q_r[UID::VEL+0]-Q_l[UID::VEL+0], Q_r[UID::VEL+1]-Q_l[UID::VEL+1], Q_r[UID::VEL+2]-Q_l[UID::VEL+2]);
	Eigen::Matrix<double, 3, 1> n2, axis1, axis2;

	if (n1(dim) < 0)
		n1 = -1*n1;
	if (n1.norm() < 1.0e-6)
		pureHLLC = true;
	else {
		n1.normalize();
		axis1 = (n1.cross(d));
		if (axis1.norm() == 0)
			pureHLL = true;
		else {
			n2 = (n1.cross(d)).cross(n1);
			if (n2.norm() == 0)
				pureHLL = true;
			else {
				n2.normalize();
				axis2 = (n2.cross(d));
				if (axis2.norm() == 0)
					pureHLLC = true;
				else
					axis2.normalize();
			}
		}
	}
	if (!pureHLLC && !pureHLL) {
		FluidArray F1, F2;
		double alpha1 = std::abs(d.dot(n1));
		double alpha2 = std::abs(d.dot(n2));
		double beta1 = 	std::sqrt(1.0 - alpha1*alpha1);
		double beta2 = 	std::sqrt(1.0 - alpha2*alpha2);
		if (axis1.norm() != 0)
			axis1.normalize();
		else
			std::cout << "ERROR: axis1 is a zero vector." << std::endl;

		Eigen::Matrix<double, 3, 3> R1 = getRotationMatrix(axis1, alpha1, beta1);
		FluidArray Q_l_R = Q_l;
		FluidArray Q_r_R = Q_r;

		rotate(R1, Q_l_R);
		rotate(R1, Q_r_R);

		m_hll.solve(F1, Q_l_R, Q_r_R, a_l2, a_r2, gamma, dim);

		Eigen::Matrix<double, 3, 3> R1_inv = getInverseRotationMatrix(axis1, alpha1, beta1);
		rotate(R1_inv, Q_l_R);
		rotate(R1_inv, Q_r_R);
		rotate(R1_inv, F1);
		Eigen::Matrix<double, 3, 3> R2 = getRotationMatrix(axis2, alpha2, beta2);
		rotate(R2, Q_l_R);
		rotate(R2, Q_r_R);

		m_hllc.solve(F2, Q_l_R, Q_r_R, a_l2, a_r2, gamma, dim);

		Eigen::Matrix<double, 3, 3> R2_inv = getInverseRotationMatrix(axis2, alpha2, beta2);
		rotate(R2_inv, Q_l_R);
		rotate(R2_inv, Q_r_R);
		rotate(R2_inv, F2);
		for (int iu = 0; iu < UID::N; ++iu)
			F[iu] = alpha1*F1[iu] + alpha2*F2[iu];
	}
	else if (pureHLL)
		m_hll.solve(F, Q_l, Q_r, a_l2, a_r2, gamma, dim);
	else if (pureHLLC)
		m_hllc.solve(F, Q_l, Q_r, a_l2, a_r2, gamma, dim);
}

std::unique_ptr<RiemannSolver> RiemannSolverFactory::create(const std::string& type, int ndims) {
	if (type.compare("HLLC") == 0)
		return std::unique_ptr<RiemannSolver>(new HartenLaxLeerContactSolver(ndims));
	else if (type.compare("RotatedHLLC") == 0)
		return std::unique_ptr<RiemannSolver>(new RotatedHartenLaxLeerSolver(ndims));
	else if (type.compare("HLL") == 0)
		return std::unique_ptr<RiemannSolver>(new HartenLaxLeerSolver(ndims));
	else if (type.compare("default") == 0)
		return std::unique_ptr<RiemannSolver>(new RotatedHartenLaxLeerSolver(ndims));
	else
		throw std::runtime_error("RiemannSolverFactory::create: unknown type " + type);
}
