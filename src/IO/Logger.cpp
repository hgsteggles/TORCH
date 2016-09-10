#include "Logger.hpp"

#include <fstream>
#include <iostream>

void FileLogPolicy::open_ostream() {
	m_outStream->open( filename.c_str(), std::ios_base::binary | std::ios_base::out );
	if ( !m_outStream->is_open())
		throw(std::runtime_error("LOGGER: Unable to open an output stream " + filename));
}


void FileLogPolicy::close_ostream() {
	if ( m_outStream ) {
		m_outStream->close();
		m_outStream->clear();
	}
}

static std::string getTime() {
	std::string time_str;
	time_t raw_time;
	time( & raw_time );
	char buffer[20];
	struct tm* timeinfo = localtime( &raw_time );
	strftime( buffer, 20, "%F %T.", timeinfo );
	time_str = buffer;
	//without the newline character
	return time_str.substr( 0 , time_str.size() );
}

std::string FileLogPolicy::getLoglineHeader() {
	std::stringstream header;
	header.str("");
	header.fill('0');
	header.width(7);
	header << m_logLineNumber++ <<" < "<< getTime() << " - ";
	header.fill('0');
	header.width(10);
	header << clock() << " > ~ ";
	return header.str();
}

void FileLogPolicy::write(const std::string& msg) {
	(*m_outStream) << getLoglineHeader() << msg << std::flush;
}

FileLogPolicy::~FileLogPolicy() {
	if ( m_outStream )
		close_ostream();
}

ConsoleLogPolicy::~ConsoleLogPolicy() {

}

void ConsoleLogPolicy::open_ostream() {

}

void ConsoleLogPolicy::close_ostream() {

}

void ConsoleLogPolicy::write(const std::string& msg) {
	(*m_outStream) << msg << std::flush;
}

Logger::~Logger() {
	for (auto& it : m_policies)
		it.second->close_ostream();
}

void Logger::registerLogPolicy(const std::string& name, std::unique_ptr<LogPolicyInterface> policy) {
	if (policy != nullptr) {
		policy->open_ostream();
		m_policies[name] = std::move(policy);
	}
}

void Logger::unregisterLogPolicy(const std::string& name) {
	m_policies.erase(name);
}

void Logger::setLogLevel(SeverityType severity) {
	for (auto& it : m_policies)
		it.second->setLogLevel(severity);
}

void Logger::setLogLevel(SeverityType severity, const std::string& policy) {
	auto it = m_policies.find(policy);
	if (it != m_policies.end())
		it->second->setLogLevel(severity);
}

void Logger::print_impl(SeverityType severity, std::unique_ptr<LogPolicyInterface>& policy) {
	if (static_cast<int>(severity) <= policy->getLogLevel())
		policy->write(m_logStream.str() );
}



