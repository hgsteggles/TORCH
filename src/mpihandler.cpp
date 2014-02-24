#include "mpihandler.hpp"

MPIHandler::MPIHandler() {
	MPI_Init(NULL, NULL);
	int name_length;
	char cname[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name(cname, &name_length);
	std::string name(cname, name_length);
	proc_name = name;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
}
MPIHandler::~MPIHandler() {
	MPI_Finalize();
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

void MPIHandler::send(double* S, int count, int destination, int tag) {
	MPI_Send(S, count, MPI_DOUBLE, destination, tag, MPI_COMM_WORLD);
}
void MPIHandler::send(int* S, int count, int destination, int tag) {
	MPI_Send(S, count, MPI_INT, destination, tag, MPI_COMM_WORLD);
}
void MPIHandler::receive(double* R, int count, int source, int tag) {
	MPI_Status status;
	MPI_Recv(R, count, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &status);
}
void MPIHandler::receive(int* R, int count, int source, int tag) {
	MPI_Status status;
	MPI_Recv(R, count, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
}

void MPIHandler::exchange(double* B, int count, int destination, int tag) {
	if (rank%2 == 0) {
		send(B, count, destination, tag);
		receive(B, count, destination, tag);
	}
	else {
		double* C = new double[count];
		memcpy( (void*)C, (void*)B, count * sizeof(double) );
		receive(B, count, destination, tag);
		send(C, count, destination, tag);
		delete[] C;
	}
}

double MPIHandler::minimum(double& x) {
	double result;
	MPI_Allreduce(&x, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return result; 
}

void MPIHandler::barrier() {
	MPI_Barrier(MPI_COMM_WORLD);
}


/* BOOST VERSION
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
*/
 */

