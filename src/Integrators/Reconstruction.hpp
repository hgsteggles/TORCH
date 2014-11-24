/*
 * reconstruction.h
 *
 *  Created on: 3 Nov 2014
 *      Author: harry
 */

#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_

class ReconstructionMethod {
public:
	virtual ~ReconstructionMethod() { }
	virtual void reconstruct(FluidArray& Q_l, FluidArray& Q_c, FluidArray& Q_r, FluidArray& left_interp, FluidArray& right_interp) const = 0;
};

class PiecewiseConstant : public ReconstructionMethod {
public:
	virtual void reconstruct(FluidArray& Q_l, FluidArray& Q_c, FluidArray& Q_r, FluidArray& left_interp, FluidArray& right_interp) const;
};

class PiecewiseLinear : public ReconstructionMethod {
public:
	virtual void reconstruct(FluidArray& Q_l, FluidArray& Q_c, FluidArray& Q_r, FluidArray& left_interp, FluidArray& right_interp) const;
	void setSlopeLimiter(std::unique_ptr<SlopeLimiter> slopeLimiter);
private:
	std::unique_ptr<SlopeLimiter> m_slopeLimiter = SlopeLimiterFactory::create("falle");
};

void PiecewiseConstant::reconstruct(FluidArray& Q_l, FluidArray& Q_c, FluidArray& Q_r, FluidArray& left_interp, FluidArray& right_interp) const {
	for (int iq = 0; iq < UID::N; ++iq) {
		left_interp[iq] = Q_c[iq];
		right_interp[iq] = Q_c[iq];
	}
}

void PiecewiseLinear::reconstruct(FluidArray& Q_l, FluidArray& Q_c, FluidArray& Q_r, FluidArray& left_interp, FluidArray& right_interp) const {
	for (int iq = 0; iq < UID::N; ++iq) {
		double dQdr = m_slopeLimiter->calculate(Q_c[iq] - Q_l[iq], Q_r[iq] - Q_c[iq]);
		left_interp[iq] = Q_c[iq] - 0.5*dQdr;
		right_interp[iq] = Q_c[iq] + 0.5*dQdr;
	}
}


class ReconstructionMethodFactory {
public:
	static std::unique_ptr<ReconstructionMethod> create(std::string type);
};

std::unique_ptr<ReconstructionMethod> ReconstructionMethodFactory::create(std::string type) {
	if (type.compare("plm") == 0)
		return std::unique_ptr<SlopeLimiter>(new PiecewiseConstant());
	else if (type.compare("pcm") == 0)
		return std::unique_ptr<SlopeLimiter>(new PiecewiseLinear());
	else {
		std::cout << "ReconstructionFactory::create: unknown type, creating PiecewiseLinear.\n";
		return std::unique_ptr<SlopeLimiter>(new PiecewiseLinear());
	}
}



#endif /* RECONSTRUCTION_H_ */
