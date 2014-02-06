/**
 * Provides the MPIHandler class.
 * @file mpihandler.hpp
 *
 * @author Harrison Steggles
 *
 * @date 28/01/2014, the first version.
 */

#ifndef MPIHANDLER_HPP_
#define MPIHANDLER_HPP_

#include <boost/mpi.hpp>
#include <boost/archive/tmpdir.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>

/**
 * @class MPIHandler
 * @brief A handler for MPI communication.
 * @version 0.3, 29/01/2014
 */
class MPIHandler {
public:
	unsigned int rank; //!< Rank of the process.
	unsigned int nproc; //!< Number of processors running this program.
	std::string proc_name; //!< Name of the processor.
	boost::mpi::communicator world; //!< MPI communicator.
	boost::mpi::environment env; //!< MPI environment.
	std::vector<boost::mpi::request> reqs; //!< Vector of MPI receive requests.

	/**
	 * @brief MPIHandler constructor.
	 * Initializes the MPI environment.
	 */
	MPIHandler();
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
	template<typename T>
	void send(int destination, int tag, T& S) {
		world.isend(destination, tag, S);
	}
	/**
	 * @brief Sends array packet to another processor.
	 * @param destination Rank of the receiving processor.
	 * @param tag Message identification tag.
	 * @param S Data to be sent.
	 * @param n Number of elements in array.
	 */
	template<typename T>
	void send(int destination, int tag, T* S, int n) {
		world.isend(destination, tag, S, n);
	}
	/**
	 * @brief Receives packet from another processor.
	 * @param source Rank of the sending processor.
	 * @param tag Message identification tag.
	 * @param R Object to receive data.
	 */
	template<typename T>
	void receive(int source, int tag, T& R){
		reqs.push_back(world.irecv(source, tag, R));
	}
	/**
	 * @brief Receives array packet from another processor.
	 * @param source Rank of the sending processor.
	 * @param tag Message identification tag.
	 * @param R Object to receive data.
	 * @param n Number of elements in array.
	 */
	template<typename T>
	void receive(int source, int tag, T* R, int n){
		reqs.push_back(world.irecv(source, tag, R, n));
	}
	/**
	 * @brief Broadcasts object from a processor to all others.
	 * @param value Object to send (if sender) or receive (if receiver).
	 * @param fromRank Rank of sending processor.
	 */
	template<typename T>
	void broadcast(T& value, int fromRank) {
		boost::mpi::broadcast(world, value, fromRank);
	}
	/**
	 * @brief Gathers objects from all other processors.
	 * @param value Object to send.
	 * @param vec Vector in which gathered objects are stored.
	 * @param toRank Rank of the receiving processor.
	 */
	template<typename T>
	void gather(T& value, std::vector<T>& vec, int toRank) {
		boost::mpi::gather(world, value, vec, toRank);
	}
	/**
	 * @brief Changes x to the minimum x passed in by all processors.
	 * @param x A value to be changed to minimum across all processors.
	 */
	void minimum(double& x);
	/**
	 * @brief Forces processor to wait until all receive requests have been fulfilled.
	 */
	void wait();
	/**
	 * @brief Forces processors to stop here until all processors reach this point.
	 */
	void barrier();
};


#endif /* MPIHANDLER_HPP_ */
