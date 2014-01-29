#include "mpihandler.hpp"

MPIHandler::MPIHandler() {
	proc_name =  env.processor_name();
	rank = world.rank();
	nproc = world.size();
}
void MPIHandler::wait() {
	if (reqs.size() > 0) {
		boost::mpi::wait_all(reqs.begin(), reqs.end());
		reqs.clear();
	}
}
int MPIHandler::getRank() {
	return rank;
}
int MPIHandler::nProcessors() {
	return nproc;
}
std::string MPIHandler::pname() {
	return proc_name;
}
std::string MPIHandler::cname() {
	std::ostringstream os;
	os << proc_name << " (" << rank << "): ";
	return os.str();
}

void MPIHandler::minimum(double& x) {
	double minval = x;
	if (rank == 0) {
		boost::mpi::reduce(world, x, minval, boost::mpi::minimum<double>(), 0);
		x = minval;
	}
	else
		boost::mpi::reduce(world, x, boost::mpi::minimum<double>(), 0);
	boost::mpi::broadcast(world, x, 0);
}

void MPIHandler::barrier() {
	world.barrier();
}

