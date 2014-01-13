#ifndef GRID3D_H
#define GRID3D_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "parameters.hpp"
#include "gridcell.hpp"

using namespace std;

class Grid3D
{
public:
	static int s_total;
	GridCell *fcell, *lcell, *fsrc;
	Boundary *fboundary;
	GridJoin *fjoin[3];
	double dx[3];
	Geometry GEOMETRY;
	int ND, NU, NR, NCELLS[3], ORDER_S, ORDER_T;
	vector<vector<vector<GridCell*> > > gridcells;

	Grid3D();
	Grid3D(const GridParameters& gp);
	~Grid3D();

	void link(int dim, GridCell* lcptr, GridCell* rcptr);
	void addSource(int x, int y, int z);
	void addGridCell(GridCell* newcptr);
	void addGhostCells(int dim, int nogcells, GridCell* cptr);
	void addBoundaryToList(Boundary* bptr);
	void addJoinToList(GridJoin* jptr, int dim);
	void addGridCellToList(GridCell* cptr);
	Boundary* createBoundary(int nogcells, int dim);
	void deleteBoundary(Boundary* bptr);

	GridCell* start();
	GridCell* traverse1D(int dim, int dc, GridCell* cptr);
	GridCell* traverse3D(int d1, int d2, int d3, int dc1, int dc2, int dc3, GridCell* cptr);
	GridCell* locate(int x, int y, int z);
	GridCell* nextCell2D(int plane, GridCell* cptr);
	GridCell* nextCell3D(GridCell* cptr);
	GridJoin* nextGridlink(int dim, GridJoin* jptr);

	GridCell* next(GridCell* srcptr, GridCell* cptr, int dxc, int dyc);
	GridCell* causalNext(GridCell* cptr, GridCell* srcptr);
	GridCell* nextSnake(GridCell* srcptr, GridCell* cptr, int dxc, int dyc, int dyz);
	GridCell* nextCausal(GridCell* cptr, GridCell* srcptr);
	double area_join(double xj[], int dim);
	void vol_cell(GridCell* cptr);

	void dice();
	void ListOpenCells();
	void ListOpenBoundaries();
};

#endif
