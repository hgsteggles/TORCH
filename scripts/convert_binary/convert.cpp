/**
 * @file convert.cpp
 *
 *  @author "Harrison Steggles"
 *  @date
 */

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cout << "Usage: " << argv[0] << " <IN_FILENAME> <OUT_FILENAME=\"conv_\"+IN_FILENAME>\n";
		return 0;
	}
	double* buffer;
	buffer = new double;
	std::ifstream fin(argv[1], std::ios_base::in | std::ios_base::binary);
	std::ostringstream oss;
	//oss << argv[1] <<
	std::string filename(argv[1]);
	std::string outfilename;
	if (argc > 2)
		outfilename = std::string(argv[2]);
	else
		outfilename = "converted_" + filename;
	std::ofstream fout(outfilename.c_str(), std::ios_base::out | std::ios_base::binary);
	int *code, *type, *ncols, *ndata;
	code = new int;
	type = new int;
	ncols = new int;
	ndata = new int;
	fin.read((char*)code, sizeof(int));
	if (*code == 211289)
		std::cout << "Acceptable Format." << std::endl;
	fin.read((char*)type, sizeof(int));
	fin.read((char*)ncols, sizeof(int));
	fin.read((char*)ndata, sizeof(int));
	fout << std::setprecision(10) << std::fixed;
	while (!fin.eof()){
		for (int i = 0; i < (*ncols)-1; ++i) {
			fin.read((char*)buffer, sizeof(double));
			if (!fin.eof())
				fout << *buffer << "\t";
		}
		fin.read((char*)buffer, sizeof(double));
		if (!fin.eof())
			fout << *buffer << "\n";
	}
	delete buffer;
	fin.close();
	fout.close();
	return 0;
}
