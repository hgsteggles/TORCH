/**
 * @file integrator.cpp
 */

#include "integrator.hpp"
#include "grid3d.hpp"
#include "radiation.hpp"
#include "hydro.hpp"
#include "io.hpp"
#include "star.hpp"
#include "gridcell.hpp"
#include "boundary.hpp"
#include "parameters.hpp"

#include <chrono>
#include <cmath>

Integrator::Integrator(IntegrationParameters& ipar, GridParameters& gpar, RadiationParameters& rpar,
		HydroParameters& hpar, PrintParameters& ppar, Scalings& scalings, MPIHandler& mpih) :
		iparams(ipar), mpihandler(mpih) {
	grid = new Grid3D(iparams.ORDER_S, gpar, scalings, mpih);
	rad = new Radiation(rpar, scalings, grid);
	hydro = new HydroDynamics(hpar, grid);
	io = new InputOutput(ppar, scalings);

	steps = 0;
}
Integrator::~Integrator() {
	delete io;
	delete hydro;
	delete rad;
	delete grid;
}

void Integrator::init() {
	steps = 0;
	Star star(0, 0, 0, 13.6*EV2ERG/grid->scale->E);
	bool snapToFace[6];
	snapToFace[0] = false;
	snapToFace[1] = false;
	snapToFace[2] = false;
	snapToFace[3] = false;
	snapToFace[4] = false;
	snapToFace[5] = false;
	rad->addStar(star, mpihandler, snapToFace);
	for(GridCell* cptr = grid->fcell; cptr != NULL; cptr = cptr->next){
		////////////
		/** STROMGREN
		double rho = rad->rparams->NHI*H_MASS_G;
		double pre = (GAS_CONST/1.0)*rho*rad->rparams->TMIN*(io->P_SCALE/io->RHO_SCALE);
		cptr->Q[ihii] = 0.0;
		cptr->Q[iden] = rho / io->RHO_SCALE;
		cptr->Q[ipre] = pre / io->P_SCALE;
		*/
		////////////
		/** SPITZER
		double rho = rad->rparams->NHI*H_MASS_G;
		double pre = (GAS_CONST/1.0)*rho*rad->rparams->THI;
		cptr->Q[ihii] = 0.0;
		double n_H = rad->rparams->NHI / (1.0/(rad->scale->L*rad->scale->L*rad->scale->L));
		double RSinf = std::pow((3.0*rad->rparams->SOURCE_S)/(4.0*PI*n_H*n_H*rad->rparams->ALPHA_B), 1.0/3.0);
		if ((cptr->xc[0]+rad->stars[0].mod[0])*grid->dx[0] <  RSinf) {
			cptr->Q[ihii] = 1.0;
			pre = (GAS_CONST/0.5)*rho*rad->rparams->THII;
		}
		cptr->Q[iden] = rho/rad->scale->RHO;
		cptr->Q[ipre] = pre/rad->scale->P;
		*/
		////////////
		/** SOD
		if ((cptr->xc[0]-(255/2.0))*(cptr->xc[0]-(255/2.0)) + (cptr->xc[1]-(255/2.0))*(cptr->xc[1]-(255/2.0)) < 1600) {
			cptr->Q[iden] = 1.0;
			cptr->Q[ipre] = 1.0;
		}
		else {
			cptr->Q[iden] = 0.125;
			cptr->Q[ipre] = 0.1;
		}
		*/
		////////////
		/** STATIONARY CONTACT
		if (cptr->xc[0] < grid->gparams->NCELLS[0]/2) {
			cptr->Q[iden] = 1.4;
			cptr->Q[ipre] = 1.0;
		}
		else {
			cptr->Q[iden] = 1.0;
			cptr->Q[ipre] = 1.0;
		}
		*/
		////////////
		/** MOVING CONTACT
		cptr->Q[ivel+0] = 0.1;
		if (cptr->xc[0] < grid->gparams->NCELLS[0]/2) {
			cptr->Q[iden] = 1.4;
			cptr->Q[ipre] = 1.0;
		}
		else {
			cptr->Q[iden] = 1.0;
			cptr->Q[ipre] = 1.0;
		}
		*/
		////////////
		/** NORMAL SHOCK
		double M = 6;
		double g = hydro->hparams->GAMMA;
		double d = 0.5;
		double F = 1.0/((2.0/(M*M*(g+1)))+((g-1)/(g+1)));
		double G = (2.0*g*M*M/((g+1)))-((g-1)/(g+1));
		double rhol = 1;
		double ul = 1;
		double pl = 1.0/(g*M*M);
		double rhor = F;
		double ur = 1.0/F;
		double pr = G/(g*M*M);
		if (cptr->xc[0] < 12) {
			cptr->Q[iden] = rhol;
			cptr->Q[ivel+0] = ul;
			cptr->Q[ipre] = pl;
		}
		else if (cptr->xc[0] == 12) {
			cptr->Q[iden] = d*rhol + (1-d)*rhor;
			cptr->Q[ivel+0] = d*ul + (1-d)*ur;
			cptr->Q[ipre] = d*pl + (1-d)*pr;
		}
		else {
			cptr->Q[iden] = rhor;
			cptr->Q[ivel+0] = ur;
			cptr->Q[ipre] = pr;
		}
		*/
		////////////
		/** DTEST */
		double rho = rad->rparams->NHI*H_MASS_G;
		double pre = (GAS_CONST/1.0)*rho*rad->rparams->THI;
		cptr->Q[ihii] = 0.0;
		double r = 0;
		for (int i = 0; i < grid->gparams->ND; ++i)
			r += (cptr->xc[i]+0.5)*(cptr->xc[i]+0.5)*grid->dx[i]*grid->dx[i];
		r = std::sqrt(r);
		double n_H = rad->rparams->NHI / (1.0/(rad->scale->L*rad->scale->L*rad->scale->L));
		//double t_rec = (1.0/(n_H*rad->rparams->ALPHA_B));
		double S = rad->rparams->SOURCE_S;
		double RSinf = pow((3.0*S)/(4.0*PI*n_H*n_H*rad->rparams->ALPHA_B), 1.0/3.0);
		if (r < RSinf) {
			cptr->Q[ihii] = 1.0;
			pre = (GAS_CONST/0.5)*rho*rad->rparams->THII;
		}
		cptr->Q[iden] = rho / io->RHO_SCALE;
		cptr->Q[ipre] = pre / io->P_SCALE;
		////////////
		/** VISHNIAC
		double rho = rad->rparams->NHI*H_MASS_G;
		double pre = (GAS_CONST/1.0)*rho*rad->rparams->THI;
		double amp = 2;
		int m = 4;
		cptr->Q[ihii] = 0.0;
		if (cptr->xc[0] < 20 + amp*std::sin(cptr->xc[1]*m*3.14159/(2.0*grid->gparams->NCELLS[1]))) {
			cptr->Q[ihii] = 1.0;
			pre = (GAS_CONST/0.5)*rho*rad->rparams->THII;
		}
		cptr->Q[iden] = rho / io->RHO_SCALE;
		cptr->Q[ipre] = pre / io->P_SCALE;
		*/
		////////////
		hydro->UfromQ(cptr->U, cptr->Q);
	}
	if(rad->nstars > 0) {
		rad->rayTrace();
		//io->printWeights(grid);
		//io->printCellPathLength(grid);
		//exit(EXIT_FAILURE);
	}
	else {
		if (mpihandler.getRank() == 0)
			std::cout << "NO SOURCES" << '\n';
	}
}

/**
 * @brief Solves the state of radiation and hydrodynamic variables until the time tmax.
 * Marches the conserved variables, Grid3D::U, in a Grid3D object until unitless (simulation) time tmax has been reached.
 * @param tmax The time (unitless) that march stops marching.
 */
void Integrator::march(const double& tmax) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed;
	double dt = 0, t = 0;
	start = std::chrono::system_clock::now();
	io->print2D(0, t, dt, grid, mpihandler);
	io->printIF(t, rad, grid, mpihandler);
	io->initProgressBar("Marching solution", mpihandler);
	while (grid->currentTime < tmax) {
		//io->freqOutput((int)(100.0*grid->currentTime/tmax), 5, rad, grid, mpihandler, "MARCH: marching solution...  ");
		io->progressBar((int)(100.0*grid->currentTime/tmax), 5, mpihandler);
		io->freqPrint(rad, grid, mpihandler);
		dt = fluidStepSplit();
		grid->currentTime += dt;
		++steps;
		//std::cerr << "t = " << grid->currentTime << ", tmax = " << tmax << ", dt = " <<  dt << '\n';
	}
	io->endProgressBar(mpihandler);
	io->print2D(100, t, dt, grid, mpihandler);
	end = std::chrono::system_clock::now();
	elapsed = end-start;
	mpihandler.barrier();
	if (mpihandler.getRank() == 0)
		std::cout << "MARCH: Took " << elapsed.count() << " seconds." << '\n';
}

/**
 * @brief Solves the state of radiation and hydrodynamic variables nsteps times.
 * Marches the conserved variables, Grid3D::U, in a Grid3D object nsteps times.
 * @param nsteps The number of steps to march.
 */
void Integrator::march(const int& nsteps) {
	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed;
	double dt = 0, t = 0;
	start = std::chrono::system_clock::now();
	io->print2D(0, t, dt, grid, mpihandler);
	io->printIF(t, rad, grid, mpihandler);
	InputOutput::initProgressBar("Marching solution", mpihandler);
	for(int step = 0; step < nsteps; grid->currentTime += dt * io->T_SCALE, ++step) {
		InputOutput::progressBar((int)(100*step/(double)nsteps), 5, mpihandler);
		io->freqPrint(rad, grid, mpihandler);
		dt = fluidStepSplit();
		++steps;
	}
	io->print2D(100, t, dt, grid, mpihandler);
	std::cout << mpihandler.cname() << "MARCH: Marching solution... " << 100 << "% " << dt << '\n';
	end = std::chrono::system_clock::now();
	elapsed = end-start;
	mpihandler.barrier();
	if (mpihandler.getRank() == 0)
		std::cout << "MARCH: Took " << elapsed.count() << " seconds." << '\n';
}
double Integrator::fluidStep() {
	double dt_hydro, dt_rad, dt, dth;
	hydro->globalQfromU();
	hydro->updateBoundaries();
	if(iparams.ORDER_S == 1)
		hydro->reconstruct();
	hydro->calcFluxes(iparams.ORDER_S);
	dt_hydro = hydro->CFL(iparams.DT_MAX);
	dt_rad = rad->calcTimeStep(iparams.DT_MAX);
	dt = std::min(dt_hydro, dt_rad);
	if (mpihandler.nProcessors() > 1)
		dt = mpihandler.minimum(dt);
	io->reduceToPrint(grid->currentTime, dt);
	grid->deltatime = dt;
	if (iparams.ORDER_T == 2) {
		dth = dt/2.0;
		hydro->globalWfromU();
	}
	else
		dth = dt;
	rad->transferRadiation(dth, mpihandler);
	//mpihandler.barrier();///////////////
	//hydro->applySrcTerms(dth);
	//resetSrcTerms();
	rad->applySrcTerms(dt, hydro->hparams);
	hydro->applySrcTerms(dth);
	hydro->advSolution(dth);
	hydro->fixSolution();
	if (iparams.ORDER_T == 2) {
		hydro->globalQfromU();
		hydro->updateBoundaries();
		//mpihandler.barrier();///////////////
		if(iparams.ORDER_S == 1)
			hydro->reconstruct();
		hydro->calcFluxes(iparams.ORDER_S);
		rad->transferRadiation(dt, mpihandler);
		hydro->globalUfromW();
		//mpihandler.barrier();///////////////
		//hydro->applySrcTerms(dt);
		//resetSrcTerms();
		rad->applySrcTerms(dt, hydro->hparams);
		hydro->applySrcTerms(dt);
		hydro->advSolution(dt);
		hydro->fixSolution();
	}
	return dt;
}

double Integrator::calcTimeStep() {
	double dt_hydro = hydro->CFL(iparams.DT_MAX);
	double dt_rad = rad->calcTimeStep(iparams.DT_MAX);
	double dt = std::min(dt_hydro, dt_rad);
	if (mpihandler.nProcessors() > 1)
		dt = mpihandler.minimum(dt);
	io->reduceToPrint(grid->currentTime, dt);
	grid->deltatime = dt;
	return dt;
}

double Integrator::fluidStepSplitOld() {
	calcHydroFlux();
	double dt = calcTimeStep();
	if(iparams.ORDER_T == 1) {
		rad->transferRadiation(dt, mpihandler);
		rad->updateSrcTerms(dt, hydro->hparams);
		advSolution(dt);
		fixSolution();
		resetSrcTerms();
		hydro->updateSrcTerms(dt);
		advSolution(dt);
		fixSolution();
		resetSrcTerms();
	}
	else {
		hydro->globalWfromU();
		rad->transferRadiation(dt/2.0, mpihandler);
		rad->updateSrcTerms(dt/2.0, hydro->hparams);
		advSolution(dt/2.0);
		fixSolution();
		resetSrcTerms();
		hydro->updateSrcTerms(dt/2.0);
		advSolution(dt/2.0);
		fixSolution();
		resetSrcTerms();
		/*
		rad->applySrcTerms(dt/2.0, hydro->hparams);
		hydro->applySrcTerms(dt/2.0);
		hydro->advSolution(dt/2.0);
		hydro->fixSolution();
		*/
		calcHydroFlux();
		rad->transferRadiation(dt, mpihandler);
		hydro->globalUfromW();
		rad->updateSrcTerms(dt, hydro->hparams);
		advSolution(dt);
		fixSolution();
		resetSrcTerms();
		hydro->updateSrcTerms(dt);
		advSolution(dt);
		fixSolution();
		resetSrcTerms();
		/*
		rad->applySrcTerms(dt/2.0, hydro->hparams);
		hydro->applySrcTerms(dt/2.0);
		hydro->advSolution(dt/2.0);
		hydro->fixSolution();
		*/
	}
	return dt;
}

double Integrator::fluidStepSplit() {
	calcHydroFlux();
	double dt = calcTimeStep();

	if (steps%2 == 0) {
		hydro->updateSrcTerms(dt/2.0);
		advSolution(dt/2.0);
		fixSolution();
		resetSrcTerms();

		hydro->globalQfromU();
		rad->transferRadiation(dt, mpihandler);
		rad->updateSrcTerms(dt, hydro->hparams);
		advSolution(dt);
		fixSolution();
		resetSrcTerms();

		calcHydroFlux();
		hydro->updateSrcTerms(dt/2.0);
		advSolution(dt/2.0);
		fixSolution();
		resetSrcTerms();
	}
	else {
		rad->transferRadiation(dt/2.0, mpihandler);
		rad->updateSrcTerms(dt/2.0, hydro->hparams);
		advSolution(dt/2.0);
		fixSolution();
		resetSrcTerms();

		calcHydroFlux();
		hydro->updateSrcTerms(dt);
		advSolution(dt);
		fixSolution();
		resetSrcTerms();

		hydro->globalQfromU();
		rad->transferRadiation(dt/2.0, mpihandler);
		rad->updateSrcTerms(dt/2.0, hydro->hparams);
		advSolution(dt/2.0);
		fixSolution();
		resetSrcTerms();
	}


	return dt;
}

void Integrator::calcHydroFlux() {
	hydro->globalQfromU();
	hydro->updateBoundaries();
	if(iparams.ORDER_S == 1)
		hydro->reconstruct();
	hydro->calcFluxes(iparams.ORDER_S);
}

void Integrator::fluidStepSplitHydro(const double& dt) {
	hydro->updateSrcTerms(dt);
}

void Integrator::fluidStepSplitRadiation(const double& dt) {
	rad->transferRadiation(dt, mpihandler);
	rad->updateSrcTerms(dt, hydro->hparams);
}

void Integrator::advSolution(const double& dt) {
	for (GridCell* cptr = grid->fcell; cptr != NULL; cptr = cptr->next) {
		for (int i = 0; i < NU; ++i) {
			cptr->U[i] += dt*cptr->UDOT[i];
		}
	}
}

void Integrator::fixSolution() {
	for (GridCell* cptr = grid->fcell; cptr != NULL; cptr = cptr->next) {
		cptr->U[iden] = std::max(cptr->U[iden], hydro->hparams->DFLOOR);
		double ke = 0.0;
		for(int dim = 0; dim < grid->gparams->ND; ++dim)
			ke += 0.5*cptr->U[ivel+dim]*cptr->U[ivel+dim]/cptr->U[iden];
		cptr->U[ipre] = std::max(cptr->U[ipre], hydro->hparams->PFLOOR/(hydro->hparams->GAMMA - 1.0) + ke);
		cptr->U[ihii] = std::max(std::min(cptr->U[ihii], cptr->U[iden]), 0.0);
	}
}

void Integrator::resetSrcTerms() {
	for (GridCell* cptr = grid->fcell; cptr != NULL; cptr = cptr->next) {
		for (int i = 0; i < NU; ++i) {
			cptr->UDOT[i] = 0;
		}
	}
}
