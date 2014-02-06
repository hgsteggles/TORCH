/* integrator.cpp */

#include "integrator.hpp"
#include "grid3d.hpp"
#include "rtmodule.hpp"
#include "hydro.hpp"
#include "io.hpp"
#include "star.hpp"
#include "gridcell.hpp"
#include "boundary.hpp"

#include <chrono>

Integrator::Integrator(GridParameters& gpar, RadiationParameters& rpar,
		HydroParameters& hpar, PrintParameters& ppar, Scalings& scalings, MPIHandler& mpih) : mpihandler(mpih) {
	grid = new Grid3D(gpar, mpih);
	rad = new Radiation(rpar, grid);
	hydro = new HydroDynamics(hpar, grid);
	io = new InputOutput(ppar, scalings);
}
Integrator::~Integrator() {
	delete io;
	delete hydro;
	delete rad;
	delete grid;
}

void Integrator::init() {
	Star star(0, 0, 0);
	rad->addStar(star, mpihandler);
	for(GridCell* cptr = grid->fcell; cptr != NULL; cptr = cptr->next){
		/** SPITZER
		double rho = rad->NHI*H_MASS_G;
		double pre = (GAS_CONST/1.0)*rho*rad->TMIN*(scale.P/scale.RHO);
		cptr->Q[ihii] = 0.0;
		if ((cptr->xc[0]-1)*grid->dx[0]*scale.L*CM2PC <  1.24989) {
			cptr->Q[ihii] = 1.0;
			pre = (GAS_CONST/0.5)*rho*rad->TMAX*(scale.P/scale.RHO);
		}
		cptr->Q[iden] = rho/scale.RHO;
		cptr->Q[ipre] = pre/scale.P;
		if(cptr == grid->fsrc)
			cptr->Q[ihii] = 1.0;
		////////////
		*/
		/** SOD
		if (cptr->xc[0] < grid->TOTNCELLS[0]/2) {
			cptr->Q[iden] = 1.0;
			cptr->Q[ipre] = 1.0;
		}
		else {
			cptr->Q[iden] = 0.125;
			cptr->Q[ipre] = 0.1;
		}
		////////////
		*/
		/** DTEST */
		double rho = rad->NHI*H_MASS_G;
		double pre = (GAS_CONST/1.0)*rho*rad->TMIN*(io->P_SCALE/io->RHO_SCALE);
		cptr->Q[ihii] = 0.0;
		double r = 0;
		for (int i = 0; i < grid->ND; ++i)
			r += (cptr->xc[i] + 0.5)*(cptr->xc[i] + 0.5)*grid->dx[i]*grid->dx[i];
		r = sqrt(r);
		if (r < (0.314*PC2CM / io->L_SCALE)) {
			cptr->Q[ihii] = 1.0;
			pre = (GAS_CONST/0.5)*rho*rad->TMAX*(io->P_SCALE/io->RHO_SCALE);
		}
		cptr->Q[iden] = rho / io->RHO_SCALE;
		cptr->Q[ipre] = pre / io->P_SCALE;
		////////////
		hydro->UfromQ(cptr->U, cptr->Q);
	}
	if(rad->nstars > 0)
		rad->rayTrace();
	else
		std::cout << "NO SOURCES" << '\n';
}

void Integrator::march(const double& tmax) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed;
	double dt = 0, t = 0;
	start = std::chrono::system_clock::now();
	std::cout << mpihandler.cname() << "MARCH: Marching solution...   0%" << '\n';
	io->print2D(0, t, dt, grid, mpihandler);
	io->printIF(t, rad, grid, mpihandler);
	io->initProgressBar("Marching solution");
	for(int prc = 0; grid->currentTime < tmax; grid->currentTime += dt){
		//io->freqOutput((int)(100.0*grid->currentTime/tmax), 5, rad, grid, mpihandler, "MARCH: marching solution...  ");
		io->progressBar((int)(100.0*grid->currentTime/tmax), 5, mpihandler);
		io->freqPrint(rad, grid, mpihandler);
		dt = fluidStep();
	}
	io->progressBar(100, 5, mpihandler);
	io->print2D(100, t, dt, grid, mpihandler);
	end = std::chrono::system_clock::now();
	elapsed = end-start;
	mpihandler.barrier();
	std::cout << "MARCH: Took " << elapsed.count() << " seconds." << '\n';
}
void Integrator::march(const int& nsteps) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed;
	double dt = 0, t = 0;
	start = std::chrono::system_clock::now();
	std::cout << mpihandler.cname() << "MARCH: Marching solution...   0%" << '\n';
	io->print2D(0, t, dt, grid, mpihandler);
	io->printIF(t, rad, grid, mpihandler);
	io->initProgressBar("Marching solution");
	for(int step = 0; step < nsteps; grid->currentTime += dt * io->T_SCALE, ++step){
		io->progressBar(step/(double)nsteps, 5, mpihandler);
		io->freqPrint(rad, grid, mpihandler);
		dt = fluidStep();
	}
	io->print2D(100, t, dt, grid, mpihandler);
	std::cout << mpihandler.cname() << "MARCH: Marching solution... " << 100 << "% " << dt << '\n';
	end = std::chrono::system_clock::now();
	elapsed = end-start;
	mpihandler.barrier();
	std::cout << "MARCH: Took " << elapsed.count() << " seconds." << '\n';
}
double Integrator::fluidStep(){
	double dt_hydro, dt_rad, dt, dth, IF = 0;
	hydro->globalQfromU();
	hydro->updateBoundaries();
	if(grid->ORDER_S == 1)
		hydro->reconstruct();
	hydro->calcFluxes();
	dt_hydro = hydro->CFL();
	dt_rad = rad->getTimeStep(dt_hydro);
	dt = std::min(dt_hydro, dt_rad);
	if (mpihandler.nProcessors() > 1)
		mpihandler.minimum(dt);
	io->reduceToPrint(grid->currentTime, dt);
	grid->deltatime = dt;
	if(grid->ORDER_T == 2){
		dth = dt/2.0;
		hydro->globalWfromU();
	}
	else
		dth = dt;
	if (mpihandler.nProcessors() == 1)
		rad->transferRadiation(dth, IF);
	else
		rad->transferRadiation(dth, IF, mpihandler);
	hydro->applySrcTerms(dth, rad);
	hydro->advSolution(dth);
	hydro->fixSolution();
	if(grid->ORDER_T == 2){
		hydro->globalQfromU();
		hydro->updateBoundaries();
		if(grid->ORDER_S == 1)
			hydro->reconstruct();
		hydro->calcFluxes();
		hydro->globalUfromW();
		if (mpihandler.nProcessors() == 1)
			rad->transferRadiation(dt, IF);
		else
			rad->transferRadiation(dt, IF, mpihandler);
		hydro->applySrcTerms(dt, rad);
		hydro->advSolution(dt);
		hydro->fixSolution();
	}
	return dt;
}
