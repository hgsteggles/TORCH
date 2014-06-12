/**
 * @file mpihandler.cpp
 */

#include "mpihandler.hpp"
#include <sstream>

/**
 * @brief MPIHandler constructor.
 * Initializes the MPI environment.
 */
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

/**
 * @brief Gets rank of the processor.
 * @return rank.
 */
int MPIHandler::getRank() {
	return rank;
}

/**
 * @brief Gets number of processors running this program.
 * @return number of processors
 */
int MPIHandler::nProcessors() {
	return nproc;
}

/**
 * @brief Gets name of processor.
 * @return name of processor.
 */
std::string MPIHandler::pname() {
	return proc_name;
}

/**
 * @brief Gets name and rank of processor.
 * @return Processor name and rank "NAME (RANK): ".
 */
std::string MPIHandler::cname() {
	std::ostringstream os;
	os << proc_name << " (" << rank << "): ";
	return os.str();
}

/**
 * @brief Sends packet to another processor.
 * @param S Message to send.
 * @param count Number of data.
 * @param destination Rank of receiving processor.
 * @param tag Message identification tag.
 */
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

void MPIHandler::write(char* filename, void* inputbuffer, int ncols, int nrows, int buffsize, BuffType btype) {
	MPI_Datatype mpitype = MPI_INTEGER;
	int typesize = 0;
	if (btype == INTEGER) {
		mpitype = MPI_INT;
		typesize = sizeof(int);
	}
	else if (btype == FLOAT) {
		mpitype = MPI_FLOAT;
		typesize = sizeof(float);
	}
	else if (btype == DOUBLE) {
		mpitype = MPI_DOUBLE;
		typesize = sizeof(double);
	}
	MPI_File thefile;
	MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
	std::string datarep = "native";
	if (rank == 0) {
		MPI_File_set_view(thefile, 0, MPI_INT, MPI_INT, (char*)datarep.c_str(), MPI_INFO_NULL);
		int buff[4] = {211289, (int)btype, ncols, nrows};
		MPI_File_write(thefile, buff, 4, MPI_INT, MPI_STATUS_IGNORE);
	}
	else
		MPI_File_set_view(thefile, (4*sizeof(int))+(rank*buffsize*typesize), MPI_INT, MPI_INT, (char*)datarep.c_str(), MPI_INFO_NULL);
	MPI_File_write(thefile, inputbuffer, buffsize, mpitype, MPI_STATUS_IGNORE);
	//barrier();
	MPI_File_close(&thefile);
}

/**
 * @brief Changes x to the minimum x passed in by all processors.
 * @param x A value to be changed to minimum across all processors.
 */
double MPIHandler::minimum(double& x) {
	double result;
	MPI_Allreduce(&x, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return result; 
}

/**
 * @brief Forces processors to stop here until all processors reach this point.
 */
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

