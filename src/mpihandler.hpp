/**
 * Provides the MPIHandler class.
 * @file mpihandler.hpp
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
 * @brief A handler for MPI communication.
 * @version 0.5, 21/02/2014
 */
class MPIHandler {
public:
	int rank; //!< Rank of the process.
	int nproc; //!< Number of processors running this program.
	std::string proc_name; //!< Name of the processor.
	

	/**
	 * @brief MPIHandler constructor.
	 * Initializes the MPI environment.
	 */
	MPIHandler();
	~MPIHandler();
	/**
	 * @brief Gets rank of the processor.
	 * @return rank.
	 */
	int getRank();
	/**
	 * @brief Gets number of processors running this program.
	 * @return number of processors
	 */
	int nProcessors();
	/**
	 * @brief Gets name of processor.
	 * @return name of processor.
	 */
	std::string pname();
	/**
	 * @brief Returns string including processor name and rank.
	 * @return Processor name and rank.
	 */
	std::string cname();
	/**
	 * @brief Sends packet to another processor.
	 * @param destination Rank of the receiving processor.
	 * @param tag Message identification tag.
	 * @param S Data to be sent.
	 */
	void send(double* S, int count, int destination, int tag);
	void send(int* S, int count, int destination, int tag);
	void receive(double* R, int count, int source, int tag);
	void receive(int* R, int count, int source, int tag);
	void exchange(double* B, int count, int destination, int tag);
	/**
	 * @brief Changes x to the minimum x passed in by all processors.
	 * @param x A value to be changed to minimum across all processors.
	 */
	double minimum(double& x);
	/**
	 * @brief Forces processors to stop here until all processors reach this point.
	 */
	void barrier();
};


#endif /* MPIHANDLER_HPP_ */
