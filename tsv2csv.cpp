#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

int main(int argc, char* argv[]) {
	std::vector<std::vector<double> > input;
	std::string line;
	std::ifstream ifile("data2D_40.txt");
	while (std::getline(ifile, line)) {
		std::vector<double> vLine;
		std::istringstream iss(line);
		double val;
		while (iss >> val)
			vLine.push_back(val);
		if (vLine.size() != 0)
			input.push_back(vLine);
	}
	ifile.close();

	std::ofstream ofile("data2D_40.csv");
	for (int i = 0; i < input.size(); ++i) {
		if (input.size() != 0)
			ofile << input[i][0];
		for (int j = 1; j < input[i].size(); ++j) {
			ofile << "," << input[i][j];
		}
		ofile << '\n';
	}
	ofile.close();
 	return 0;
}
