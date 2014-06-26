/**
 * @file io.cpp
 */

#include "io.h"
#include "gridcell.h"
#include "star.h"

#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/bzip2.hpp>

/**
 * @brief InputOutput constructor.
 * @param rp Parameters to pass in.
 * @param sc Scalings to pass in for printing physical units.
 */
InputOutput::InputOutput(PrintParameters& pp, Scalings& sc) : scale(sc) {
	DIR_2D = pp.DIR_2D;
	DIR_IF = pp.DIR_IF;
	PRINT2D_ON = pp.PRINT2D_ON;
	PRINTIF_ON = pp.PRINTIF_ON;
}

int InputOutput::percent = 0;
bool InputOutput::freqPrinting = false;
std::string InputOutput::progressMessage = "Progress";
bool InputOutput::debugging = false;

void InputOutput::initProgressBar(std::string msg, MPIHandler& mpihandler) {
	InputOutput::percent = 0;
	InputOutput::freqPrinting = false;
	InputOutput::progressMessage = msg;
}
void InputOutput::progressBar(const int& prcnow, const int& every, MPIHandler& mpihandler) {
	if (prcnow - InputOutput::percent >= every) {
		if (!InputOutput::debugging && mpihandler.getRank() == 0) {
			int barWidth = 40;
			int pos = barWidth * prcnow / 100.0;
			std::cout << InputOutput::progressMessage << ": [";
			for (int i = 0; i < barWidth; ++i) {
				if (i < pos) std::cout << "=";
				else if (i == pos) std::cout << ">";
				else std::cout << " ";
			}
			std::cout << "] " << prcnow << " %" << '\r';
			std::cout << std::flush;
		}
		InputOutput::percent = prcnow;
		InputOutput::freqPrinting = true;
		if (InputOutput::percent == 100 && mpihandler.getRank() == 0)
			std::cout << std::endl;
	}
	else
		InputOutput::freqPrinting = false;
}
void InputOutput::endProgressBar(MPIHandler& mpihandler){
	InputOutput::progressBar(100, 1, mpihandler);
	InputOutput::percent = 0;
	InputOutput::freqPrinting = false;
	InputOutput::progressMessage = "Progress";
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
	if (mpih.getRank() != 0) {
		int ready;
		mpih.receive(&ready, 1, mpih.getRank() - 1, PRINTSTARBENCH_MSG);
	}
	{
		std::ostringstream os;
		os << DIR_2D << "starbench_" << i << ".txt.bz2";
		std::ofstream ofile(os.str().c_str(), std::ios_base::app | std::ios_base::binary);
		if(!ofile)
			std::cerr << "ERROR: Unable to open " << os.str() << '\n';
		boost::iostreams::filtering_ostream out;
		out.push(boost::iostreams::bzip2_compressor());
		out.push(ofile);
		out << std::setprecision(10) << std::fixed;
		for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
			double x1 = 0, x2 = 0, v1 = 0, v2 = 0;
			if (gptr->gparams->ND > 1) {
				x1 = scale.fromCodeUnits((cptr->xc[1]+0.5)*gptr->dx[1], 0, 1, 0)*CM2PC;
				v1 = scale.fromCodeUnits(cptr->Q[ivel+1], 0, 1, -1)*0.001;
			}
			else if (gptr->gparams->ND > 2) {
				x2 = scale.fromCodeUnits((cptr->xc[2]+0.5)*gptr->dx[2], 0, 1, 0)*CM2PC;
				v2 = scale.fromCodeUnits(cptr->Q[ivel+2], 0, 1, -1)*0.001;
			}
			out << fortranformat(scale.fromCodeUnits((cptr->xc[0]+0.5)*gptr->dx[0], 0, 1, 0)*CM2PC, 16, 7, 3);
			out << fortranformat(x1, 16, 7, 3);
			out << fortranformat(x2, 16, 7, 3);
			out << fortranformat(scale.fromCodeUnits(cptr->Q[ivel+0], 0, 1, -1)*0.001, 16, 7, 3);
			out << fortranformat(v1, 16, 7, 3);
			out << fortranformat(v2, 16, 7, 3);
			out	<< fortranformat(scale.fromCodeUnits(cptr->Q[iden], 1, -3, 0), 16, 7, 3);
			out << fortranformat((rad->rparams.THI + cptr->Q[ihii]*(rad->rparams.THII - rad->rparams.THI)), 16, 7, 3);
			out << fortranformat(cptr->Q[ihii], 16, 7, 3);
			out << '\n';
		}
	}
	if (mpih.getRank() != mpih.nProcessors() - 1) {
		int ready = 1;
		mpih.send(&ready, 1, mpih.getRank() + 1, PRINTSTARBENCH_MSG);
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

void InputOutput::printBinary2D(const int& step, const double& t, const double& dt, Grid3D* gptr, MPIHandler& mpih) const {
	if (PRINT2D_ON) {
		std::ostringstream os;
		os << DIR_2D << "data2D_";
		os << step << ".hbin";
		int ncols = 5;
		if (gptr->gparams->ND > 1)
			ncols += 2;
		if (gptr->gparams->ND > 2)
			ncols += 2;
		const int nrows = gptr->gparams->NCELLS[0]*gptr->gparams->NCELLS[1]*gptr->gparams->NCELLS[2];
		const int nbuff = gptr->gparams->CORECELLS[0]*gptr->gparams->CORECELLS[1]*gptr->gparams->CORECELLS[2]*ncols;
		double* buff = new double[nbuff];
		int i = 0;
		for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next) {
			for (int id = 0; id < gptr->gparams->ND; ++id)
				buff[i++] = cptr->xc[id];
			buff[i++] = cptr->Q[iden];
			buff[i++] = cptr->Q[ipre];
			buff[i++] = cptr->Q[ihii];
			for (int id = 0; id < gptr->gparams->ND; ++id)
				buff[i++] = cptr->Q[ivel+id];
		}
		if (i != nbuff) {
			std::cout << "ERROR: buffer not filled in printBinary2D.\n";
			exit(EXIT_FAILURE);
		}
		mpih.write((char*)os.str().c_str(), buff, ncols, nrows, nbuff, mpih.DOUBLE);
		delete[] buff;
	}
}

/**
 * @brief Prints primitive variables for a 2D slice of the simulation grid in the Grid3D object pointed to by gptr.s
 * @param step Number to append to the output filename.
 * @param t Current simulation time.
 * @param dt Current delta-time.
 * @param gptr Pointer to Grid3D object to print.
 * @param mpih Provides MPI information for printing with multiple cores.
 */
void InputOutput::print2D(const int& step, const double& t, const double& dt, Grid3D* gptr, MPIHandler& mpih) const {
	if (PRINT2D_ON) {
		if (mpih.getRank() != 0) {
			int ready;
			mpih.receive(&ready, 1, mpih.getRank() - 1, PRINT2D_MSG);
		}
		{
			/* creating filename */
			std::ostringstream os;
			os << DIR_2D << "data2D_";
			os << step << ".txt.bz2";
			/* opening new file for appending data */
			std::ofstream ofile(os.str().c_str(), std::ios_base::app | std::ios_base::binary);
			if(!ofile)
				std::cerr << "ERROR: Unable to open " << os.str() << '\n';
			boost::iostreams::filtering_ostream buff;
			buff.push(boost::iostreams::bzip2_compressor());
			buff.push(ofile);
			/* writing data to file */
			std::stringstream data;
			buff << std::setprecision(10) << std::fixed;
			for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next){
				buff << cptr->xc[0];
				if(gptr->gparams->ND > 1)
					buff << '\t' << cptr->xc[1];
				if(gptr->gparams->ND > 2)
					buff << '\t' << cptr->xc[2];
				buff << '\t' << cptr->Q[iden];
				//buff << '\t' << std::sqrt(cptr->Q[iden]*cptr->Q[ivel+0]*cptr->Q[ivel+0]/(1.4*cptr->Q[ipre]));
				buff << '\t' << cptr->Q[ipre];
				buff << '\t' << cptr->Q[ihii];
				buff << '\t' << cptr->Q[ivel+0];
				if(gptr->gparams->ND > 1)
					buff << '\t' << cptr->Q[ivel+1];
				if(gptr->gparams->ND > 2)
					buff << '\t' << cptr->Q[ivel+2];
				buff << '\n';
			}
		}
		if (mpih.getRank() != mpih.nProcessors() - 1) {
			int ready = 1;
			mpih.send(&ready, 1, mpih.getRank() + 1, PRINT2D_MSG);
		}
	}
}

/**
 * @brief Calculates and prints the location of the ionization front from a Star that is located at the origin.
 * @param t Current simulation time.
 * @param rad Radiation provides this method with radiation parameters.
 * @param gptr The grid on which the Star is located.
 * @param mpih Provides MPI information for calculating and printing with multiple cores.
 */
void InputOutput::printIF(const double& t, const Radiation* rad, Grid3D* gptr, MPIHandler& mpih) const {
	if (PRINTIF_ON && rad->nstars > 0) {
		/* IONIZATION FRONT DETECTION */
		double IF = 0;
		bool found = false;
		if (mpih.getRank() != 0) {
			int found_int;
			mpih.receive(&found_int, 1, mpih.getRank() - 1, PRINTIF_NEXT_MSG);
			found = found_int == 1 ? true : false;
		}
		if(!found) {
			for(GridCell* cptr = gptr->fcell; cptr != NULL; cptr = gptr->traverse1D(0, 1, cptr)) {
				double frac = cptr->U[ihii]/cptr->U[iden];
				if(cptr->xc[0] >= gptr->fcell->xc[0] && frac < 0.5 && !found) { // IF at HII fraction = 0.5 .
					if(cptr != gptr->fcell){
						if(cptr->left[0] != NULL) {
							double fracLeft = cptr->left[0]->U[ihii]/cptr->left[0]->U[iden];
							double interp = (0.5-fracLeft)/(frac-fracLeft);
							IF = (cptr->left[0]->xc[0] + rad->stars[0].mod[0]) - rad->stars[0].x[0] + interp;
						}
						else
							IF = (cptr->xc[0] + rad->stars[0].mod[0]) - rad->stars[0].x[0];
						//IF = sqrt(IF*IF + 0.25);
					}
					found = true;
					break;
				}
			}
			if(mpih.getRank() != mpih.nProcessors()-1) {
				int found_int = found ? 1 : 0;
				mpih.send(&found_int, 1, mpih.getRank()+1, PRINTIF_NEXT_MSG);
			}
			if(mpih.getRank() != 0) {
				int found_int = found ? 1 : 0;
				mpih.send(&found_int, 1, 0, PRINTIF_FOUND_MSG);
				if(found)
					mpih.send(&IF, 1, 0, PRINTIF_IF_MSG);
			}
			else if(!found) {
				for(int i = 1; i < mpih.nProcessors(); i++) {
					int found_int = found ? 1 : 0;
					mpih.receive(&found_int, 1, i, PRINTIF_FOUND_MSG);
					if(found) {
						mpih.receive(&IF, 1, i, PRINTIF_IF_MSG);
						break;
					}
				}
			}
		}
		if(mpih.getRank() == 0) {
			IF /= gptr->gparams->NCELLS[0];
			std::ofstream outFile("tmp/IF.dat", std::ios::app);
			double n_H, S = 0, RSinf;
			double cII, RI;
			n_H = scale.toCodeUnits(rad->rparams.NHI, 0, -3, 0);
			//t_rec = (1.0/(n_H*rad->rparams.ALPHA_B));
			if (rad->stars.size() != 0)
				S = rad->stars[0].photonRate;
			RSinf = pow((3.0*S)/(4.0*PI*n_H*n_H*rad->rparams.ALPHA_B), 1.0/3.0);
			cII = scale.toCodeUnits(sqrt(GAS_CONST*2.0*rad->rparams.THII), 0, 1, -1);
			//RI = pow((1.0-exp(-t/t_rec)), 1.0/3.0)*RSinf;
			RI = RSinf*pow(1+(7*cII*t)/(4*RSinf), 4.0/7.0);
			double t_s = RSinf/cII;
			double error = 0;
			if (RI != 0.0)
				error = fabs(IF-RI)/RI;
			outFile << t/t_s << '\t';
			outFile << scale.fromCodeUnits(IF, 0, 1, 0)*CM2PC << '\t';
			outFile << scale.fromCodeUnits(RI, 0, 1, 0)*CM2PC << '\t';
			outFile << error << '\n';
			outFile.close();
		}
	}
}
void InputOutput::printWeights(Grid3D* gptr) {
	std::ofstream ofile("tmp/weights.dat", std::ios::app);
	for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next) {
		ofile << "{ ";
		for (int i = 0; i < 3; ++i) {
			if (cptr->xc[i] < 10)
				ofile << " " << cptr->xc[i] << " ";
			else
				ofile << cptr->xc[i] << " ";
		}
		ofile << "}:   weights of { ";

		ofile << std::fixed << std::setprecision(3);
		for (int i = 0; i < 4; ++i) {
			std::ostringstream os;
			if (cptr->NN[i] != NULL) {
				double diff = cptr->NN[i]->xc[0]-cptr->xc[0];
				if (diff > 0)
					os << "(+x)";
				else if (diff < 0)
					os << "(-x)";
				diff = cptr->NN[i]->xc[1]-cptr->xc[1];
				if (diff > 0)
					os << "(+y)";
				else if (diff < 0)
					os << "(-y)";
				diff = cptr->NN[i]->xc[2]-cptr->xc[2];
				if (diff > 0)
					os << "(+z)";
				else if (diff < 0)
					os << "(-z)";
			}
			else {
				os << "(##)";
				if (i > 0)
					os << "(##)";
				if (i > 2)
					os << "(##)";
			}
			ofile << os.str() << " ";
		}
		ofile << "} = { ";
		for (int i = 0; i < 4; ++i)
			ofile << cptr->NN_weights[i] << " ";
		ofile << "}\n";
	}
	ofile.close();
}

void InputOutput::printCellPathLength(Grid3D* gptr) {
	std::ofstream ofile("tmp/ds.dat", std::ios::app);
	for (GridCell* cptr = gptr->fcell; cptr != NULL; cptr = cptr->next) {
		ofile << "{ ";
		for (int i = 0; i < 3; ++i) {
			if (cptr->xc[i] < 10)
				ofile << " " << cptr->xc[i] << " ";
			else
				ofile << cptr->xc[i] << " ";
		}
		ofile << "}:   ds = ";
		ofile << std::fixed << std::setprecision(5);
		ofile << cptr->ds << '\n';
	}
	ofile.close();
}

/**
 * @brief Prints all parameters associated with this simulation.
 * TO-DO: Write function that provides consistent parameter output.
 * @param runtime Simulation time.
 * @param gpar Grid3D parameters for printing.
 * @param rpar Radiation parameters for printing.
 * @param hpar Hydrodynamics parameters for printing.
 * @param ppar I/O parameters for printing.
 * @param scale Scales for printing.
 */
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
