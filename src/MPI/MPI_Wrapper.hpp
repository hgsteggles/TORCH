/** Provides the MPIW class.
 *
 * @file MPI_Wrapper.hpp
 *
 * @author Harrison Steggles
 */

#ifndef MPIHANDLER_HPP_
#define MPIHANDLER_HPP_

#include <functional>
#include <string>

enum class SendID : unsigned int {PARTITION_MSG, RADIATION_MSG, THERMO_MSG, PRINT2D_MSG,
	CFL_COLLECT, CFL_BROADCAST, PRINTIF_NEXT_MSG, PRINTIF_FOUND_MSG,
	PRINTIF_IF_MSG, PRINTSTARBENCH_MSG, PRINT_HEATING_MSG, PERIODIC_MSG};
enum BuffType {INTEGER, FLOAT, DOUBLE}; //!< buffer data types.

/**
 * @class MPIW
 *
 * @brief A singleton class for MPI communication.
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
	int rank = 0; //!< Rank of the process.
	int nproc = 1; //!< Number of processors running this program.
	std::string proc_name = ""; //!< Name of the processor.

	~MPIW();

	// Getters/Setters.
	int getRank() const;
	int nProcessors() const;
	std::string pname() const;
	std::string cname() const;

	// Message passing methods.
	void send(double* S, int count, int destination, SendID tag) const;
	void send(int* S, int count, int destination, SendID tag) const;
	void receive(double* R, int count, int source, SendID tag) const;
	void receive(int* R, int count, int source, SendID tag) const;
	void exchange(double* send_buffer, double* receive_buffer, int count, int destination, SendID tag) const;
	void barrier() const;
	void broadcastBoolean(bool msg, int source) const;
	void broadcastString(std::string& msg, int source) const;

	// Output.
	void write(char* filename, void* inputbuffer, int ncols, int nrows, int buffsize, BuffType btype) const;

	// Misc. methods.
	double minimum(double& x) const;
	double maximum(double& x) const;
	double sum(double& x) const;

	void serial(const std::function<void()>& f);
	void abort();

private:
    MPIW(int* argc, char*** argv);
    MPIW(MPIW const&);
    void operator=(MPIW const&);
};


#endif // MPIHANDLER_HPP_
