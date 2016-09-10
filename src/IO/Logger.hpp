/** Provides the Logger class.
 *
 * @file Logger.hpp
 *
 * @author Harrison Steggles
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
#include <iostream>
#include <map>

enum class SeverityType : unsigned int { FATAL_ERROR = 1, ERROR, WARNING, NOTICE, INFO, DEBUG };

/**
 * @class LogPolicyInterface
 * @brief Interface skeleton for any logging policy.
 */
class LogPolicyInterface {
public:
	virtual ~LogPolicyInterface() {};
	virtual void open_ostream() = 0;
	virtual void close_ostream() = 0;
	virtual void write(const std::string& msg) = 0;
	virtual int getLogLevel() final {return logLevel;};
	virtual void setLogLevel(SeverityType level) final {logLevel = static_cast<int>(level);};
private:
	int logLevel = static_cast<int>(SeverityType::DEBUG);
};

/**
 * @class FileLogPolicy
 * @brief A logger with this policy will output logs to a file.
 */
class FileLogPolicy : public LogPolicyInterface {
public:
	FileLogPolicy(const std::string& filename) : filename(filename), m_outStream( new std::ofstream ) {};
	~FileLogPolicy();
	virtual void open_ostream();
	virtual void close_ostream();
	virtual void write(const std::string& msg);

	std::string getLoglineHeader();
private:
	unsigned m_logLineNumber;
	std::string filename = "";
	std::unique_ptr< std::ofstream > m_outStream;
};

/**
 * @class ConsoleLogPolicy
 * @brief A logger with this policy will output logs to standard output.
 */
class ConsoleLogPolicy : public LogPolicyInterface {
public:
	ConsoleLogPolicy() : m_outStream( &std::cout ) {};
	~ConsoleLogPolicy();
	void open_ostream();
	void close_ostream();
	void write(const std::string& msg);
private:
	std::ostream* m_outStream;
};

/**
 * @class Logger
 *
 * @brief Handles all logging information.
 */
class Logger {
public:
	static Logger& Instance() {
		static Logger instance;
		return instance;
	}

	void setLogLevel(SeverityType severity);
	void setLogLevel(SeverityType severity, const std::string& policy);
	void registerLogPolicy(const std::string& name, std::unique_ptr<LogPolicyInterface> policy);
	void unregisterLogPolicy(const std::string& name);

	~Logger();

	template< SeverityType severity , typename...Args >
	void print( Args...args );
private:
	SeverityType logLevel = SeverityType::DEBUG;
	std::stringstream m_logStream;
	std::map<std::string, std::unique_ptr<LogPolicyInterface>> m_policies;
	std::mutex m_writeMutex;

	Logger( ) {};
	Logger( Logger const& );
	void operator=( Logger const& );

	//Core printing functionality
	void print_impl(SeverityType severity, std::unique_ptr<LogPolicyInterface>& policy);
	template<typename First, typename...Rest>
	void print_impl(SeverityType severity, std::unique_ptr<LogPolicyInterface>& policy, First parm1, Rest...parm);
};

template< SeverityType severity , typename...Args >
void Logger::print( Args...args ) {
	m_writeMutex.lock();

	switch( severity ) {
		case SeverityType::DEBUG:
			m_logStream << "<DEBUG>: ";
			break;
		case SeverityType::INFO:
			m_logStream << "<INFO>: ";
			break;
		case SeverityType::NOTICE:
			m_logStream << "<NOTICE>: ";
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
	for (auto& it : m_policies) {
		print_impl(severity, it.second, args...);
	}
	m_logStream.str("");

	m_writeMutex.unlock();
}

template<typename First, typename...Rest >
void Logger::print_impl(SeverityType severity, std::unique_ptr<LogPolicyInterface>& policy, First parm1, Rest...parm) {
	m_logStream << parm1;
	print_impl(severity, policy, parm...);
}



#endif // LOGGER_HPP_
