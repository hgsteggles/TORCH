#include "MPI_Wrapper.hpp"

#include <stdlib.h>
#include <iostream>
#include <sstream>

#include <mpi.h>

/**
 * @brief MPIHandler constructor.
 * Initializes the MPI environment.
 */
MPIW::MPIW(int* argc, char*** argv) {
	MPI_Init(argc, argv);
	int name_length = 0;
	char cname[MPI_MAX_PROCESSOR_NAME];
	MPI_Get_processor_name(cname, &name_length);
	std::string name(cname, name_length);
	proc_name = name;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
}
MPIW::~MPIW() {
	MPI_Finalize();
}

/**
 * @brief Gets rank of the processor.
 * @return rank.
 */
int MPIW::getRank() const {
	return rank;
}

/**
 * @brief Gets number of processors running this program.
 * @return number of processors
 */
int MPIW::nProcessors() const {
	return nproc;
}

/**
 * @brief Gets name of processor.
 * @return name of processor.
 */
std::string MPIW::pname() const {
	return proc_name;
}

/**
 * @brief Gets name and rank of processor.
 * @return Processor name and rank "NAME (RANK): ".
 */
std::string MPIW::cname() const {
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
void MPIW::send(double* S, int count, int destination, SendID tag) const {
	MPI_Send((void*)S, count, MPI_DOUBLE, destination, (int)tag, MPI_COMM_WORLD);
}
void MPIW::send(int* S, int count, int destination, SendID tag) const {
	MPI_Send(S, count, MPI_INT, destination, (int)tag, MPI_COMM_WORLD);
}
void MPIW::receive(double* R, int count, int source, SendID tag) const {
	MPI_Status status;
	MPI_Recv((void*)R, count, MPI_DOUBLE, source, (int)tag, MPI_COMM_WORLD, &status);
}
void MPIW::receive(int* R, int count, int source, SendID tag) const {
	MPI_Status status;
	MPI_Recv(R, count, MPI_INT, source, (int)tag, MPI_COMM_WORLD, &status);
}

void MPIW::exchange(double* send_buffer, double* receive_buffer, int count, int destination, SendID tag) const {
	if (rank%2 == 0 || (nproc == 2 && tag == SendID::PERIODIC_MSG && rank == 1)) {
		send(send_buffer, count, destination, tag);
		receive(receive_buffer, count, destination, tag);
	}
	else {
		receive(receive_buffer, count, destination, tag);
		send(send_buffer, count, destination, tag);
	}
}

void MPIW::write(char* filename, void* inputbuffer, int ncols, int nrows, int buffsize, BuffType btype) const {
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
double MPIW::minimum(double& x) const {
	double result;
	MPI_Allreduce(&x, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	return result; 
}

/**
 * @brief Changes x to the maximum x passed in by all processors.
 * @param x A value to be changed to maximum across all processors.
 */
double MPIW::maximum(double& x) const {
	double result;
	MPI_Allreduce(&x, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	return result;
}

double MPIW::sum(double& x) const {
	double result;
	MPI_Allreduce(&x, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return result;
}

/**
 * @brief Forces processors to stop here until all processors reach this point.
 */
void MPIW::barrier() const {
	MPI_Barrier(MPI_COMM_WORLD);
}

void MPIW::broadcastBoolean(bool msg, int source) const {
	MPI_Bcast(&msg, 1, MPI_INT, source, MPI_COMM_WORLD);
}

void MPIW::broadcastString(std::string& msg, int source) const {
	int length;
	if (rank == source)
		length = msg.size();
	MPI_Bcast(&length, 1, MPI_INT, source, MPI_COMM_WORLD);
	char* cmsg;
	if (rank == source)
		cmsg = (char*) msg.c_str();
	else
		cmsg = (char*) malloc(sizeof(char)*(length + 1));
	MPI_Bcast(cmsg, length, MPI_CHAR, source, MPI_COMM_WORLD);
	if (rank != source)
		msg = std::string(cmsg);
}

void MPIW::serial(const std::function<void()>& f) {
	for (int iproc = 0; iproc < nproc; ++iproc) {
		if (rank == iproc)
			f();
		barrier();
	}
}

void MPIW::abort() {
	MPI_Abort(MPI_COMM_WORLD, 0);
}
