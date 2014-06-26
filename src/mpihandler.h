/**
 * Provides the MPIHandler class.
 * @file mpihandler.h
 *
 * @author Harrison Steggles
 *
 * @date 28/01/2014, the first version.
 * @date 21/02/2014 - mpi implementation overhauled. Uses mpi functions directly rather than uses boost interface.
 * This solved memory problems when sending large packets of data between core boundaries.
 */

#ifndef MPIHANDLER_HPP_
#define MPIHANDLER_HPP_

#include "mpi.h"

/**
 * @class MPIHandler
 *
 * @brief A handler for MPI communication.
 *
 * @version 0.7, 13/06/2014
 */
class MPIHandler {
public:
	int rank; //!< Rank of the process.
	int nproc; //!< Number of processors running this program.
	std::string proc_name; //!< Name of the processor.

	enum BuffType {INTEGER, FLOAT, DOUBLE}; //!< buffer data types.

	//Structors.
	MPIHandler();
	~MPIHandler();

	//Setters and Getters.
	int getRank();
	int nProcessors();
	std::string pname();
	std::string cname();

	//Message passing methods.
	void send(double* S, int count, int destination, int tag);
	void send(int* S, int count, int destination, int tag);
	void receive(double* R, int count, int source, int tag);
	void receive(int* R, int count, int source, int tag);
	void exchange(double* B, int count, int destination, int tag);
	void barrier();

	//Output.
	void write(char* filename, void* inputbuffer, int ncols, int nrows, int buffsize, BuffType btype);

	//Misc. methods.
	double minimum(double& x);
};


#endif /* MPIHANDLER_HPP_ */
