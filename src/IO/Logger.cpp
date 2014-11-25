#include "Logger.hpp"

void FileLogPolicy::open_ostream(const std::string& name) {
	m_outStream->open( name.c_str(), std::ios_base::binary | std::ios_base::out );
	if ( !m_outStream->is_open() )
		throw(std::runtime_error("LOGGER: Unable to open an output stream"));
}

void FileLogPolicy::close_ostream() {
	if ( m_outStream )
		m_outStream->close();
}

void FileLogPolicy::write(const std::string& msg) {
	(*m_outStream) << msg << std::endl;
}

FileLogPolicy::~FileLogPolicy() {
	if ( m_outStream )
		close_ostream();
}


