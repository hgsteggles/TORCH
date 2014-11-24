/**
 * @file io.cpp
 */

#include "DataPrinter.hpp"
#include "Fluid.hpp"
//#include "grid3d.h"
#include "GridFactory.hpp"
#include "Radiation.hpp"
#include "Hydro.hpp"
#include "Constants.hpp"
#include "GridCell.hpp"
#include "Star.hpp"
#include "Converter.hpp"

#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string.h>

#include "boost/iostreams/filtering_stream.hpp"
#include "boost/iostreams/filtering_streambuf.hpp"
#include "boost/iostreams/filter/bzip2.hpp"

void print(std::string str) {
	std::cout << str << std::endl;
}

void DataPrinter::initialise(std::shared_ptr<Constants> c, std::string output_directory) {
	consts = std::move(c);
	dir2D = output_directory;
	printing_on = dir2D != "";
}

void DataPrinter::printSTARBENCH(const Radiation& rad, const Hydrodynamics& hydro, const Fluid& fluid) {
	MPIW& mpihandler = MPIW::Instance();
	const Grid& grid = fluid.getGrid();
	for(int i = (int)printTimes.size()-1; i >= 0 ; --i) {
		if (std::abs(grid.currentTime - printTimes[i])/grid.currentTime <= 0.0000000001 && !printDone[i]) {
			fluid.globalQfromU();
			if (mpihandler.getRank() != 0) {
				int ready;
				mpihandler.receive(&ready, 1, mpihandler.getRank() - 1, SendID::PRINTSTARBENCH_MSG);
			}
			{
				std::ostringstream os;
				os << dir2D << "starbench_" << i << ".txt.bz2";
				std::ofstream ofile(os.str().c_str(), std::ios_base::app | std::ios_base::binary);
				if(!ofile)
					std::cerr << "ERROR: Unable to open " << os.str() << '\n';
				boost::iostreams::filtering_ostream out;
				out.push(boost::iostreams::bzip2_compressor());
				out.push(ofile);
				out << std::setprecision(10) << std::fixed;
				for (GridCell& cell : grid.getCells()){
					double x1 = 0, x2 = 0, v1 = 0, v2 = 0;
					if (consts->nd > 1) {
						x1 = consts->converter.CM_2_PC(consts->converter.fromCodeUnits(cell.xc[1]*grid.dx[1], 0, 1, 0));
						v1 = consts->converter.fromCodeUnits(cell.Q[UID::VEL+1], 0, 1, -1)*0.001;
					}
					if (consts->nd > 2) {
						x2 = consts->converter.CM_2_PC(consts->converter.fromCodeUnits(cell.xc[2]*grid.dx[2], 0, 1, 0));
						v2 = consts->converter.fromCodeUnits(cell.Q[UID::VEL+2], 0, 1, -1)*0.001;
					}
					out << fortranformat(consts->converter.CM_2_PC(consts->converter.fromCodeUnits(cell.xc[0]*grid.dx[0], 0, 1, 0)), 16, 7, 3);
					out << fortranformat(x1, 16, 7, 3);
					out << fortranformat(x2, 16, 7, 3);
					out	<< fortranformat(consts->converter.fromCodeUnits(cell.Q[UID::DEN], 1, -3, 0), 16, 7, 3);
					out << fortranformat(cell.Q[UID::HII], 16, 7, 3);
					out << fortranformat(consts->converter.fromCodeUnits(cell.Q[UID::PRE], 1, -1, -2), 16, 7, 3);
					out << fortranformat((rad.THI + cell.Q[UID::HII]*(rad.THII - rad.THI)), 16, 7, 3);
					out << fortranformat(consts->converter.fromCodeUnits(cell.Q[UID::VEL+0], 0, 1, -1)*0.001, 16, 7, 3);
					out << fortranformat(v1, 16, 7, 3);
					out << fortranformat(v2, 16, 7, 3);
					out << '\n';
				}
			}
			if (mpihandler.getRank() != mpihandler.nProcessors() - 1) {
				int ready = 1;
				mpihandler.send(&ready, 1, mpihandler.getRank() + 1, SendID::PRINTSTARBENCH_MSG);
			}
			printDone[i] = true;
		}
	}
}

std::string DataPrinter::fortranformat(double value, int w, int d, int e) const {
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

/*
boost::iostreams::filtering_ostream createCompressionStream(std::string directory, int step) {
	// Creating filename.
	std::ostringstream os;
	os << directory << "data2D_";
	os << step << ".txt.bz2";

	std::ofstream ofile(os.str().c_str(), std::ios_base::app | std::ios_base::binary);
	if (!ofile)
		std::cerr << "ERROR: Unable to open " << os.str() << '\n';
	boost::iostreams::filtering_ostream buff;
	buff.push(boost::iostreams::bzip2_compressor());
	buff.push(ofile);
	return buff;
}
*/
void DataPrinter::printBinary2D(const int step, const double t, const Grid& grid) const {
	if (!printing_on)
		return;
	std::ostringstream os;
	os << dir2D << "data2D_";
	os << step << ".hbin";
	int ncols = 5;
	if (consts->nd > 1)
		ncols += 2;
	if (consts->nd > 2)
		ncols += 2;
	const int nrows = grid.ncells[0]*grid.ncells[1]*grid.ncells[2];
	const int nbuff = grid.coreCells[0]*grid.coreCells[1]*grid.coreCells[2]*ncols;
	double* buff = new double[nbuff];
	int i = 0;
	for (GridCell& cell : grid.getCells()) {
		for (int id = 0; id < consts->nd; ++id)
			buff[i++] = cell.xc[id];
		buff[i++] = cell.Q[UID::DEN];
		buff[i++] = cell.Q[UID::PRE];
		buff[i++] = cell.Q[UID::HII];
		for (int id = 0; id < consts->nd; ++id)
			buff[i++] = cell.Q[UID::VEL+id];
	}
	if (i != nbuff)
		throw std::runtime_error("DataPrinter::printBinary2D: buffer not filled in printBinary2D.");

	MPIW::Instance().write((char*)os.str().c_str(), buff, ncols, nrows, nbuff, BuffType::DOUBLE);
	delete[] buff;
}

/**
 * @brief Prints primitive variables for a 2D slice of the simulation grid in the Grid object pointed to by grid.s
 * @param step Number to append to the output filename.
 * @param t Current simulation time.
 * @param dt Current delta-time.
 * @param scale For conversion from unit-less variables to variables in cgs units.
 * @param grid The Grid to print.
 * @param mpih Provides MPI information for printing with multiple cores.
 */
void DataPrinter::print2D(const std::string& append_name, const double t, const Grid& grid) const {
	MPIW& mpihandler = MPIW::Instance();
	if (!printing_on)
		return;
	mpihandler.serial([&] () {
		// Creating filename.
		std::ostringstream os;
		os << dir2D << "data2D_";
		os << append_name << ".txt.bz2";

		std::ofstream ofile(os.str().c_str(), std::ios_base::app | std::ios_base::binary);
		if(!ofile)
			std::cerr << "ERROR: Unable to open " << os.str() << '\n';
		boost::iostreams::filtering_ostream buff;
		buff.push(boost::iostreams::bzip2_compressor());
		buff.push(ofile);

		// Writing data to file.
		std::stringstream data;
		buff << std::setprecision(10) << std::scientific;
		if (mpihandler.getRank() == 0) {
			buff << consts->converter.fromCodeUnits(t, 0, 0, 1) << '\n';
			buff << grid.ncells[0] << '\n';
			buff << grid.ncells[1] << '\n';
			buff << grid.ncells[2] << '\n';
		}
		for (GridCell& cell : grid.getCells()) {
			for (int idim = 0; idim < consts->nd; ++idim)
				buff << consts->converter.fromCodeUnits(cell.xc[idim]*grid.dx[idim], 0, 1, 0) << '\t';
			buff << consts->converter.fromCodeUnits(cell.Q[UID::DEN], 1, -3, 0);
			buff << '\t' << consts->converter.fromCodeUnits(cell.Q[UID::PRE], 1, -1, -2);
			buff << '\t' << cell.Q[UID::HII];
			for (int idim = 0; idim < consts->nd; ++idim)
				buff << '\t' << consts->converter.fromCodeUnits(cell.Q[UID::VEL+idim], 0, 1, -1);
			buff << '\n';
		}
	});
}

/**
 * @brief Calculates and prints the location of the ionization front from a Star that is located at the origin.
 * @param t Current simulation time.
 * @param scale For conversion from unit-less variables to variables in cgs units.
 * @param rad Radiation provides this method with radiation parameters.
 * @param grid The grid on which the Star is located.
 * @param mpih Provides MPI information for calculating and printing with multiple cores.
 *
 * @note TODO FIX THIS.
 */
/*
void InputOutput::printIF(const double t, const Radiation& rad, const Grid& grid, const MPIHandler& mpih) const {
	if (PRINTIF_ON && rad.stars.size() > 0 && mpihandler.nProcessors() == 0) {
		// IONIZATION FRONT DETECTION
		double IF = 0;
		bool found = false;
		if (mpihandler.getRank() != 0) {
			int found_int;
			mpihandler.receive(&found_int, 1, mpihandler.getRank() - 1, SendID::PRINTIF_NEXT_MSG);
			found = found_int == 1 ? true : false;
		}
		if(!found) {
			for(GridCell* cptr = grid.fcell; cptr != NULL; cptr = grid.traverse1D(0, 1, cptr)) {
				double frac = cell.U[ihii]/cell.U[iden];
				if(cell.xc[0] >= grid.fcell->xc[0] && frac < 0.5 && !found) { // IF at HII fraction = 0.5 .
					if(cptr != grid.fcell){
						if(cell.left[0] != NULL) {
							double fracLeft = cell.left[0]->U[ihii]/cell.left[0]->U[iden];
							double interp = (0.5-fracLeft)/(frac-fracLeft);
							IF = cell.left[0]->xc[0] - rad.stars[0].xc[0] + interp;
						}
						else
							IF = cell.xc[0] - rad.stars[0].xc[0];
						//IF = sqrt(IF*IF + 0.25);
					}
					found = true;
					break;
				}
			}
			if(mpihandler.getRank() != mpihandler.nProcessors()-1) {
				int found_int = found ? 1 : 0;
				mpihandler.send(&found_int, 1, mpihandler.getRank()+1, SendID::PRINTIF_NEXT_MSG);
			}
			if(mpihandler.getRank() != 0) {
				int found_int = found ? 1 : 0;
				mpihandler.send(&found_int, 1, 0, SendID::PRINTIF_FOUND_MSG);
				if(found)
					mpihandler.send(&IF, 1, 0, SendID::PRINTIF_IF_MSG);
			}
			else if(!found) {
				for(int i = 1; i < mpihandler.nProcessors(); i++) {
					int found_int = found ? 1 : 0;
					mpihandler.receive(&found_int, 1, i, SendID::PRINTIF_FOUND_MSG);
					if(found) {
						mpihandler.receive(&IF, 1, i, SendID::PRINTIF_IF_MSG);
						break;
					}
				}
			}
		}
		if(mpihandler.getRank() == 0) {
			IF /= grid.ncells[0];
			std::ofstream outFile("tmp/IF.dat", std::ios::app);
			double n_H, S = 0, RSinf;
			double cII, RI;
			n_H = consts->converter.toCodeUnits(rad.nHI, 0, -3, 0);
			//t_rec = (1.0/(n_H*rad.ALPHA_B));
			if (rad.stars.size() != 0)
				S = rad.stars[0].sparams.photonRate;
			RSinf = pow((3.0*S)/(4.0*PI*n_H*n_H*rad.alphaB), 1.0/3.0);
			cII = consts->converter.toCodeUnits(sqrt(GAS_CONST*2.0*rad.THII), 0, 1, -1);
			//RI = pow((1.0-exp(-t/t_rec)), 1.0/3.0)*RSinf;
			RI = RSinf*pow(1+(7*cII*t)/(4*RSinf), 4.0/7.0);
			double t_s = RSinf/cII;
			double error = 0;
			if (RI != 0.0)
				error = fabs(IF-RI)/RI;
			outFile << t/t_s << '\t';
			outFile << consts->converter.fromCodeUnits(IF, 0, 1, 0)*CM2PC << '\t';
			outFile << consts->converter.fromCodeUnits(RI, 0, 1, 0)*CM2PC << '\t';
			outFile << error << '\n';
			outFile.close();
		}
	}
}
 */

void DataPrinter::printHeating(const int step, const double t, const Grid& grid) const {
	MPIW& mpihandler = MPIW::Instance();
	mpihandler.serial([&] () {
		/* creating filename */
		std::ostringstream os;
		os << dir2D << "heating_";
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
		buff << std::setprecision(10) << std::scientific;
		if (mpihandler.getRank() == 0) {
			buff << consts->converter.fromCodeUnits(t, 0, 0, 1) << '\n';
			buff << grid.ncells[0] << '\n';
			buff << grid.ncells[1] << '\n';
			buff << grid.ncells[2] << '\n';
		}
		for (GridCell& cell : grid.getCells()) {
			for (int idim = 0; idim < consts->nd; ++idim)
				buff << cell.xc[idim] << '\t';
			buff << consts->converter.fromCodeUnits(cell.H[0], 1, -1, -3);
			for (int i = 1; i < HID::N; ++i)
				buff << '\t' << consts->converter.fromCodeUnits(cell.H[i], 1, -1, -3);
			buff << '\n';
		}
	});
}

void DataPrinter::printVariables(const int step, const double t, const Grid& grid) const {
	MPIW& mpihandler = MPIW::Instance();
	mpihandler.serial([&] () {
		/* creating filename */
		std::ostringstream os;
		os << dir2D << "heating_";
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
		buff << std::setprecision(10) << std::scientific;
		if (mpihandler.getRank() == 0) {
			buff << consts->converter.fromCodeUnits(t, 0, 0, 1) << '\n';
			buff << grid.ncells[0] << '\n';
			buff << grid.ncells[1] << '\n';
			buff << grid.ncells[2] << '\n';
		}
		for (GridCell& cell : grid.getCells()) {
			for (int idim = 0; idim < consts->nd; ++idim)
				buff << cell.xc[idim] << '\t';
			buff << cell.H[0];
			for (int i = 1; i < HID::N; ++i)
				buff << '\t' << cell.H[i];
			buff << '\n';
		}
	});
}

void DataPrinter::printVariable(const int step, const double t, const Grid& grid) const {
	MPIW& mpihandler = MPIW::Instance();
	mpihandler.serial([&] () {
		/* creating filename */
		std::ostringstream os;
		os << dir2D << "var_";
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
		buff << std::setprecision(10) << std::scientific;
		if (mpihandler.getRank() == 0) {
			buff << consts->converter.fromCodeUnits(t, 0, 0, 1) << '\n';
			buff << grid.ncells[0] << '\n';
			buff << grid.ncells[1] << '\n';
			buff << grid.ncells[2] << '\n';
		}
		for (GridCell& cell : grid.getCells()) {
			for (int idim = 0; idim < consts->nd; ++idim)
				buff << cell.xc[idim] << '\t';
			buff << cell.H[0] << '\n';
		}
	});
}

void DataPrinter::printWeights(const Grid& grid) const {
	std::ofstream ofile("tmp/weights.dat", std::ios::app);
	for (GridCell& cell : grid.getCells()) {
		ofile << "{ ";
		for (int i = 0; i < 3; ++i) {
			if (cell.xc[i] < 10)
				ofile << " " << cell.xc[i] << " ";
			else
				ofile << cell.xc[i] << " ";
		}
		ofile << "}:   weights of { ";

		ofile << std::fixed << std::setprecision(3);
		for (int i = 0; i < 4; ++i) {
			std::ostringstream os;
			if (cell.NN[i] != NULL) {
				double diff = cell.NN[i]->xc[0]-cell.xc[0];
				if (diff > 0)
					os << "(+x)";
				else if (diff < 0)
					os << "(-x)";
				diff = cell.NN[i]->xc[1]-cell.xc[1];
				if (diff > 0)
					os << "(+y)";
				else if (diff < 0)
					os << "(-y)";
				diff = cell.NN[i]->xc[2]-cell.xc[2];
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
			ofile << cell.NN_weights[i] << " ";
		ofile << "}\n";
	}
	ofile.close();
}
/*
void InputOutput::printCellPathLength(const Grid& grid) const {
	std::ofstream ofile("tmp/ds.dat", std::ios::app);
	for (GridCell* cptr = grid.fcell; cptr != NULL; cptr = cell.next) {
		ofile << "{ ";
		for (int i = 0; i < 3; ++i) {
			if (cell.xc[i] < 10)
				ofile << " " << cell.xc[i] << " ";
			else
				ofile << cell.xc[i] << " ";
		}
		ofile << "}:   ds = ";
		ofile << std::fixed << std::setprecision(5);
		ofile << cell.ds << '\n';
	}
	ofile.close();
}
 */
void DataPrinter::fileToMap(const std::string& myString, std::map<double,double> myMap) const {
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

void DataPrinter::addPrintTime(const double t) {
	printTimes.push_back(t);
	printDone.push_back(false);
}
void DataPrinter::reduceToPrint(const double currTime, double& dt) const {
	for (int i = 0; i < (int)printTimes.size(); i++) {
		if(currTime < printTimes[i] && currTime+dt > printTimes[i]) {
			dt = std::min(dt, printTimes[i]-currTime);
		}
	}
}

PrintParameters::PrintParameters(){
	dir2D = "tmp/";
}

void PrintParameters::printInfo() const {
	std::cout << "dir2D = " << dir2D << "\n";
}
