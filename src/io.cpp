/* io.cpp */

#include "io.hpp"
#include "gridcell.hpp"
#include "star.hpp"

#include <fstream>
#include <cmath>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

InputOutput::InputOutput(const PrintParameters& pp, const Scalings& sc) {
	DIR_2D = pp.DIR_2D;
	DIR_IF = pp.DIR_IF;
	PRINT2D_ON = pp.PRINT2D_ON;
	PRINTIF_ON = pp.PRINTIF_ON;
	L_SCALE = sc.L;
	M_SCALE = sc.M;
	T_SCALE = sc.T;
	V_SCALE = sc.V;
	RHO_SCALE = sc.RHO;
	P_SCALE = sc.P;
	E_SCALE = sc.E;
	percent = 0;
	freqPrinting = false;
	debugging = false;
	progressMessage = "Progress";
}
void InputOutput::initProgressBar(std::string msg) {
	percent = 0;
	freqPrinting = false;
	progressMessage = msg;
}
void InputOutput::progressBar(const int& prcnow, const int& every, MPIHandler& mpihandler) {
	if (prcnow - percent >= every) {
		if (!debugging && mpihandler.getRank() == 0) {
			int barWidth = 40;
			int pos = barWidth * prcnow / 100.0;
			std::cout << progressMessage << ": [";
			for (int i = 0; i < barWidth; ++i) {
				if (i < pos) std::cout << "=";
				else if (i == pos) std::cout << ">";
				else std::cout << " ";
			}
			std::cout << "] " << prcnow << " %" << '\r';
			std::cout << std::flush;
		}
		percent = prcnow;
		freqPrinting = true;
		if (percent == 100)
			std::cout << std::endl;
	}
	else
		freqPrinting = false;
}
void InputOutput::freqPrint(const Radiation* rad, Grid3D* gptr, MPIHandler& mpihandler) {
	if (freqPrinting) {
		if (debugging)
			std::cout << mpihandler.cname() << progressMessage << ": " << percent << "% " << gptr->deltatime << '\n';
		print2D(percent, gptr->currentTime, gptr->deltatime, gptr, mpihandler);
		printIF(gptr->currentTime, rad, gptr, mpihandler);
	}
	for(int i = 0; i < (int)printTimes.size(); ++i) {
		if (std::abs(gptr->currentTime - printTimes[i])/gptr->currentTime <= 0.0000000001)
			printSTARBENCH(i, rad, gptr, mpihandler);
	}
}
void InputOutput::printSTARBENCH(const int& i, const Radiation* rad, Grid3D* gptr, MPIHandler& mpih) {
	std::ostringstream os;
	os << DIR_2D << "starbench_" << i << ".txt.gz";
	std::ofstream ofile(os.str().c_str(), std::ios_base::app | std::ios_base::binary);
	if(!ofile)
		std::cerr << "ERROR: Unable to open " << os.str() << '\n';
	boost::iostreams::filtering_ostream out;
    out.push(boost::iostreams::gzip_compressor());
    out.push(ofile);
	out << std::setprecision(10) << std::fixed;
	//out << "t" << '\t' << "x" << '\t' << "y" << '\t' << "z" << '\t' << "HIIfrac" << '\t' << "rho" << '\t' << "pre" << '\n';
	if (mpih.getRank() != 0) {
		bool ready;
		mpih.receive(mpih.getRank() - 1, PRINTSTARBENCH_MSG, ready);
		mpih.wait();
	}
	for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
		double x1 = 0, x2 = 0, v1 = 0, v2 = 0;
		if (gptr->ND > 1) {
			x1 = (cptr->xc[1]+0.5)*gptr->dx[1]*L_SCALE*CM2PC;
			v1 = cptr->Q[ivel+1]*V_SCALE*0.001;
		}
		else if (gptr->ND > 2) {
			x2 = (cptr->xc[2]+0.5)*gptr->dx[2]*L_SCALE*CM2PC;
			v2 = cptr->Q[ivel+2]*V_SCALE*0.001;
		}
		out << fortranformat((cptr->xc[0]+0.5)*gptr->dx[0]*L_SCALE*CM2PC, 16, 7, 3);
		out << fortranformat(x1, 16, 7, 3);
		out << fortranformat(x2, 16, 7, 3);
		out << fortranformat(cptr->Q[ivel+0]*V_SCALE*0.001, 16, 7, 3);
		out << fortranformat(v1, 16, 7, 3);
		out << fortranformat(v2, 16, 7, 3);
		out	<< fortranformat(cptr->Q[iden]*RHO_SCALE, 16, 7, 3);
		out << fortranformat((rad->TMIN + cptr->Q[ihii]*(rad->TMAX - rad->TMIN))*(P_SCALE/RHO_SCALE), 16, 7, 3);
		out << fortranformat(cptr->Q[ihii], 16, 7, 3);
		out << '\n';
	}
	if (mpih.getRank() != mpih.nProcessors() - 1) {
		bool ready = true;
		mpih.send(mpih.getRank() + 1, PRINTSTARBENCH_MSG, ready);
	}
}
std::string InputOutput::fortranformat(double value, int w, int d, int e) {
	//find exponent
	double num = std::abs(value);
	int exponent = 0;
	for (; num > 1.0; exponent++) num /= 10.0;
	for (; num < 0.1 && num > 0.0; exponent--) num *= 10.0;
	// total width must be greater than e + d + 5
	if (w < e + d + 5) w = e + d + 5;
	std::ostringstream os;
	for (int i = 0; i < w - e - d - 5; ++i) os << " ";
	if (value < 0)
		os << "-";
	else
		os << " ";
	/*
	os << "0.";
	int nextnumber = (int)((num*std::pow(10, d)) + 0.5);
	if (nextnumber != 0)
		os << nextnumber;
	else {
		for (int i = 0; i < d; ++i)
			os << "0";
	}
	*/
	os << std::fixed << std::setprecision(d) << num;
	os << "E";
	if (exponent >= 0)
		os << "+";
	else
		os << "-";
	double exptmp = std::abs(exponent) + 0.5;
	int edigits = 0;
	while (exptmp > 1) {
		exptmp /= 10.0;
		++edigits;
	}
	if (edigits == 0) edigits = 1;
	for(int i = 0; i < e - edigits; ++i)
		os << "0";
	os << std::abs(exponent);
	return os.str();
}
void InputOutput::print2D(const int& step, const double& t, const double& dt, Grid3D* gptr, MPIHandler& mpih) const {
	if (PRINT2D_ON) {
		/* creating filename */
		std::ostringstream os;
		os << DIR_2D << "data2D_";
		os << step << ".txt.gz";
		/* opening new file for appending data */
		std::ofstream ofile(os.str().c_str(), std::ios_base::app | std::ios_base::binary);
		if(!ofile)
			std::cerr << "ERROR: Unable to open " << os.str() << '\n';
		boost::iostreams::filtering_ostream out;
	    out.push(boost::iostreams::gzip_compressor());
	    out.push(ofile);
		/* writing data to file */
		out << std::setprecision(10) << std::fixed;
		//out << "t" << '\t' << "x" << '\t' << "y" << '\t' << "z" << '\t' << "HIIfrac" << '\t' << "rho" << '\t' << "pre" << '\n';
		if (mpih.getRank() != 0) {
			bool ready;
			mpih.receive(mpih.getRank() - 1, PRINT2D_MSG, ready);
			mpih.wait();
		}
		for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
			out << cptr->xc[0];
			if(gptr->ND > 1)
				out << '\t' << cptr->xc[1];
			if(gptr->ND > 2)
				out << '\t' << cptr->xc[2];
			out << '\t'	<< cptr->Q[iden];
			out << '\t' << cptr->Q[ipre];
			out << '\t' << cptr->Q[ihii];
			out << '\t' << cptr->Q[ivel+0];
			if(gptr->ND > 1)
				out << '\t' << cptr->Q[ivel+1];
			if(gptr->ND > 2)
				out << '\t' << cptr->Q[ivel+2];
			out << '\n';
		}
		if (mpih.getRank() != mpih.nProcessors() - 1) {
			bool ready = true;
			mpih.send(mpih.getRank() + 1, PRINT2D_MSG, ready);
		}
	}
}
void InputOutput::printIF(const double& t, const Radiation* rad, Grid3D* gptr, MPIHandler& mpih) const {
	if (PRINTIF_ON && rad->nstars > 0) {
		/* IONIZATION FRONT DETECTION */
		double IF = 0;
		bool found = false;
		if (mpih.getRank() != 0) {
			mpih.receive(mpih.getRank() - 1, PRINTIF_NEXT_MSG, found);
			mpih.wait();
		}
		if(!found) {
			for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->traverse1D(0, 1, cptr)) {
				double frac = cptr->U[ihii]/cptr->U[iden];
				if(cptr->xc[0] >= gptr->fcell->xc[0] && frac < 0.5 && !found) { // IF at HII fraction = 0.5 .
					if(cptr != gptr->fcell){
						if(cptr->left[0] != NULL) {
							double fracLeft = cptr->left[0]->U[ihii]/cptr->left[0]->U[iden];
							double interp = (0.5-fracLeft)/(frac-fracLeft);
							IF = (cptr->left[0]->xc[0] + 0.5) - rad->stars[0].x[0] + interp;
						}
						else
							IF = (cptr->xc[0] + 0.5) - rad->stars[0].x[0];
						IF = sqrt(IF*IF + 0.25);
					}
					found = true;
					break;
				}
			}
			if(mpih.getRank() != mpih.nProcessors()-1)
				mpih.send(mpih.getRank()+1, PRINTIF_NEXT_MSG, found);
			if(mpih.getRank() != 0) {
				mpih.send(0, PRINTIF_FOUND_MSG, found);
				if(found)
					mpih.send(0, PRINTIF_IF_MSG, IF);
			}
			else if(!found) {
				for(int i = 1; i < mpih.nProcessors(); i++) {
					mpih.receive(i, PRINTIF_FOUND_MSG, found);
					mpih.wait();
					if(found) {
						mpih.receive(i, PRINTIF_IF_MSG, IF);
						mpih.wait();
						break;
					}
				}
			}
		}
		if(mpih.getRank() == 0) {
			IF /= gptr->NCELLS[0];
			std::ofstream outFile("tmp/IF.dat", std::ios::app);
			double n_H, S, RSinf;
			double cII, RI;
			//double RS, t_rec;
			n_H = rad->NHI / (1.0/(L_SCALE*L_SCALE*L_SCALE));
			//t_rec = (1.0/(n_H*ALPHA));
			double t1 = t;
			S = rad->SOURCE_S;
			RSinf = pow((3.0*S)/(4.0*PI*n_H*n_H*rad->ALPHA_B), 1.0/3.0);
			//RS = pow((1.0-exp(-t1/t_rec)), 1.0/3.0)*RSinf;
			cII = sqrt(GAS_CONST*2.0*rad->TMAX) / V_SCALE;
			RI = RSinf*pow(1+(7*cII*t1)/(4*RSinf), 4.0/7.0);
			double t_s = RSinf/cII;
			double error = 0;
			if (RI != 0.0)
				error = fabs(IF-RI)/RI;
			outFile << t/t_s << '\t';
			outFile << IF*L_SCALE*CM2PC << '\t';
			outFile << RI*L_SCALE*CM2PC << '\t';
			outFile << error << '\n';
			outFile.close();
		}
	}
}
/*
void printParams(double runtime, GridParameters& gpar, RadiationParameters& rpar, HydroParameters& hpar, PrintParameters& ppar, Scalings& scale){
	ofstream outFile("tmp/params.txt");
	outFile << setprecision(8);
	outFile << scientific;
	outFile << setw(15) << left << "dt" << '\t' << p.dt * T_SCALE << '\n';
	outFile << setw(15) << left << "nocells_x" << '\t' << p.nocells[0] << '\n';
	outFile << setw(15) << left << "nocells_y" << '\t' << p.nocells[1] << '\n';
	outFile << setw(15) << left << "nocells_z" << '\t' << p.nocells[2] << '\n';
	outFile << setw(15) << left << "x width" << '\t' << L_SCALE << '\n';
	outFile << setw(15) << left << "geom" << '\t' << p.geom << '\n';
	outFile << setw(15) << left << "scheme" << '\t' << p.scheme << '\n';
	outFile << setw(15) << left << "K1" << '\t' << K1 << '\n';
	outFile << setw(15) << left << "K2" << '\t' << K2 << '\n';
	outFile << setw(15) << left << "K3" << '\t' << K3 << '\n';
	outFile << setw(15) << left << "K4" << '\t' << K4 << '\n';
	outFile << setw(15) << left << "hv" << '\t' << p.hv * E_SCALE << '\n';
	outFile << setw(15) << left << "srcS" << '\t' << SOURCE_S * (1.0/T_SCALE) << '\n';
	outFile << setw(15) << left << "tau_0" << '\t' << TAU_0 << '\n';
	outFile << setw(15) << left << "PIcrossSect" << '\t' << PI_CROSS_SECT * (L_SCALE*L_SCALE) << '\n';
	outFile << setw(15) << left << "alphaB" << '\t' << ALPHA * (L_SCALE*L_SCALE*L_SCALE/T_SCALE) << '\n';
	outFile << setw(15) << left << "clumpBL_x" << '\t' << p.clump[0] << '\n';
	outFile << setw(15) << left << "clumpBL_y" << '\t' << p.clump[1] << '\n';
	outFile << setw(15) << left << "clumpW" << '\t' << p.clump[2] << '\n';
	outFile << setw(15) << left << "no_t_steps" << '\t' << p.tstepctr << '\n';
	outFile << setw(15) << left << "runtime" << '\t' << runtime << '\n';
	outFile << setw(15) << left << "units" << '\t' << "cgs" << '\n';
	outFile.close();
}
*/
void InputOutput::fileToMap(const std::string& myString, std::map<double,double> myMap) const {
	std::ifstream inFile(myString.c_str(), std::ios::in);
	if(!inFile)
	  std::cout << "ERROR: unable to open " << myString << '\n';
	else{
		double a,b;
		while(!inFile.eof() ){
			inFile >> a >> b;
			myMap[a] = b;
		}
	}
  inFile.close();
}
