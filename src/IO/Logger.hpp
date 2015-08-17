/** Provides the Logger class.
 *
 * @file Logger.hpp
 *
 * @author Harrison Steggles
 *
 * @date 28/01/2014 - the first version.
 */

#ifndef LOGGER_HPP_
#define LOGGER_HPP_

#include <ctime>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>

#include "../MPI/MPI_Wrapper.hpp"

/**
 * @class LogPolicyInterface
 *
 * @brief Interface skeleton for any logging policy.
 *
 * @version 0.8, 24/11/2014
 */
class LogPolicyInterface {
public:
	virtual ~LogPolicyInterface() {};
	virtual void open_ostream(const std::string& name) = 0;
	virtual void close_ostream() = 0;
	virtual void write(const std::string& msg) = 0;
};

/**
 * @class FileLogPolicy
 *
 * @brief A logger with this policy will output logs to a file.
 *
 * @version 0.8, 24/11/2014
 */
class FileLogPolicy : public LogPolicyInterface {
public:
	FileLogPolicy() : m_outStream( new std::ofstream ) {};
	~FileLogPolicy();
	void open_ostream(const std::string& name);
	void close_ostream();
	void write(const std::string& msg);
private:
	std::unique_ptr< std::ofstream > m_outStream;
};

enum class SeverityType : unsigned int { DEBUG = 1, FATAL_ERROR, ERROR, WARNING };

/**
 * @class Logger
 *
 * @brief Handles all logging information.
 *
 * @version 0.8, 24/11/2014
 */
template< typename log_policy >
class Logger {
public:
	static Logger& Instance() {
		std::string filename = "log/torch.log" + std::to_string(MPIW::Instance().getRank());
		static Logger instance( filename, true );
		return instance;
	}

	~Logger();

	template< SeverityType severity , typename...Args >
	void print( Args...args );
private:
	unsigned m_logLineNumber;
	std::string getTime();
	std::string getLoglineHeader();
	std::stringstream m_logStream;
	log_policy* m_policy;
	std::mutex m_writeMutex;

	Logger( const std::string& name, bool isRootProcess );
	Logger( Logger const& );
	void operator=( Logger const& );

	//Core printing functionality
	void print_impl();
	template<typename First, typename...Rest>
	void print_impl(First parm1, Rest...parm);
};

template< typename log_policy >
Logger< log_policy >::Logger( const std::string& name, bool isRootProcess ) {
	m_logLineNumber = 0;
	if (isRootProcess) {
		m_policy = new log_policy;
		if( !m_policy )
			throw std::runtime_error("LOGGER: Unable to create the logger instance.");
		m_policy->open_ostream( name );
	}
}

template< typename log_policy >
Logger< log_policy >::~Logger() {
	if ( m_policy ) {
		m_policy->close_ostream();
		delete m_policy;
	}
}

template< typename log_policy >
template< SeverityType severity , typename...Args >
void Logger< log_policy >::print( Args...args ) {
	if (!m_policy)
		return;
	m_writeMutex.lock();
	switch( severity )
	{
	case SeverityType::DEBUG:
		m_logStream << "<DEBUG>: ";
		break;
	case SeverityType::WARNING:
		m_logStream << "<WARNING>: ";
		break;
	case SeverityType::ERROR:
		m_logStream << "<ERROR>: ";
		break;
	case SeverityType::FATAL_ERROR:
		m_logStream << "<FATAL_ERROR>: ";
		break;
	};
	print_impl( args... );
	m_writeMutex.unlock();
}

template< typename log_policy >
void Logger< log_policy >::print_impl() {
	m_policy->write( getLoglineHeader() + m_logStream.str() );
	m_logStream.str("");
}

template< typename log_policy >
template<typename First, typename...Rest >
void Logger< log_policy >::print_impl(First parm1, Rest...parm) {
	m_logStream << parm1;
	print_impl(parm...);
}

template< typename log_policy >
std::string Logger< log_policy >::getTime() {
	std::string time_str;
	time_t raw_time;
	time( & raw_time );
	char buffer[20];
	struct tm* timeinfo = localtime( &raw_time );
	strftime( buffer, 20, "%F %T.", timeinfo );
	time_str = buffer;
	//without the newline character
	return time_str.substr( 0 , time_str.size() - 1 );
}

template< typename log_policy >
std::string Logger< log_policy >::getLoglineHeader() {
	std::stringstream header;
	header.str("");
	header.fill('0');
	header.width(7);
	header << m_logLineNumber++ <<" < "<< getTime() << " - ";
	header.fill('0');
	header.width(7);
	header << clock() << " > ~ ";
	return header.str();
}



#endif // LOGGER_HPP_
