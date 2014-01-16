/* gridcell.h */

#ifndef GRIDCELL_H
#define GRIDCELL_H

#include <stdio.h>
#include <map>
#include <vector>
#include "parameters.hpp"
using namespace std;
class Cell;
class GridJoin;
class GhostCell;
class Boundary;

class GridCell{
	public:
		GridCell();
		GridCell(int i, int j, int k);
		~GridCell();
		/* SETS */
		void set_U(int index, double value);
		void set_xcs(int x, int y, int z);
		/* GETS */
		int get_xc(int i);
		double get_U(int index);
		/* FUNCTIONS */
		double temperature();

		static int s_total;
		Boundary *bd[3];
		GridJoin *rjoin[3], *ljoin[3];
		GridCell *right[3], *left[3];
		GridCell *next;
		double U[NU];
		double Q[NU];
		double W[NU];
		double R[NR];
		int xc[3];
		double UL[3][NU];
		double UR[3][NU];
		double QL[3][NU];
		double QR[3][NU];
		double vol;
};

class GridJoin{
	public:
		GridJoin();
		~GridJoin();
		static int s_total;
		GridCell *lcell, *rcell;
		GhostCell *inghost, *outghost;
		GridJoin* next;
		double F[NU];
		double xj[NU];
		double area;
};

#endif
