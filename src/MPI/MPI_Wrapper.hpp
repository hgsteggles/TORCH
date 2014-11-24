/** Provides the MPIW class.
 *
 * @file MPI_Wrapper.hpp
 *
 * @author Harrison Steggles
 *
 * @date 28/01/2014, the first version.
 * @date 21/02/2014 - mpi implementation changed. Uses mpi functions directly rather than boost interface.
 * @date 24/11/2014 - MPIW::exchange modified to accept two buffers (instead of allocating heap space for a temporary each time this
 * is called.
 * @date 24/11/2014 - broadcastString, broadcastBoolean, serial, print, abort, sum and maximum methods added.
 */

#ifndef MPIHANDLER_HPP_
#define MPIHANDLER_HPP_

#include <mpi.h>
#include <string>
#include <functional>

enum class SendID : unsigned int {PARTITION_MSG, RADIATION_MSG, THERMO_MSG, PRINT2D_MSG,
	CFL_COLLECT, CFL_BROADCAST, PRINTIF_NEXT_MSG, PRINTIF_FOUND_MSG,
	PRINTIF_IF_MSG, PRINTSTARBENCH_MSG, PRINT_HEATING_MSG, PERIODIC_MSG};
enum BuffType {INTEGER, FLOAT, DOUBLE}; //!< buffer data types.

/**
 * @class MPIW
 *
 * @brief A handler for MPI communication.
 *
 * @version 0.8, 24/11/2014
 */
class MPIW {
public:
	static MPIW& Instance(int* argc, char*** argv) {
		static MPIW instance(argc, argv);
		return instance;
	}
	static MPIW& Instance() {
		return Instance(0, 0);
	}
	int rank; //!< Rank of the process.
	int nproc; //!< Number of processors running this program.
	std::string proc_name; //!< Name of the processor.

	~MPIW();

	//Setters and Getters.
	int getRank() const;
	int nProcessors() const;
	std::string pname() const;
	std::string cname() const;

	//Message passing methods.
	void send(double* S, int count, int destination, SendID tag) const;
	void send(int* S, int count, int destination, SendID tag) const;
	void receive(double* R, int count, int source, SendID tag) const;
	void receive(int* R, int count, int source, SendID tag) const;
	void exchange(double* send_buffer, double* receive_buffer, int count, int destination, SendID tag) const;
	void barrier() const;
	void broadcastBoolean(bool msg, int source) const;
	void broadcastString(std::string& msg, int source) const;

	//Output.
	void write(char* filename, void* inputbuffer, int ncols, int nrows, int buffsize, BuffType btype) const;

	//Misc. methods.
	double minimum(double& x) const;
	double maximum(double& x) const;
	double sum(double& x) const;

	void serial(const std::function<void()>& f);
	void print(const std::string& message);
	void abort();

private:
    MPIW(int* argc, char*** argv);
    MPIW(MPIW const&);
    void operator=(MPIW const&);
};


#endif /* MPIHANDLER_HPP_ */
