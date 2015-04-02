// ============================================================================
// gzstream, C++ iostream classes wrapping the zlib compression library.
// Copyright (C) 2001  Deepak Bandyopadhyay, Lutz Kettner
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : gzstream.h
// Revision      : $Revision: 1.5 $
// Revision_date : $Date: 2002/04/26 23:30:15 $
// Author(s)     : Deepak Bandyopadhyay, Lutz Kettner
// 
// Standard streambuf implementation following Nicolai Josuttis, "The 
// Standard C++ Library".
// ============================================================================

#ifndef STREAMGZ_HPP_
#define STREAMGZ_HPP_

#include <libio.h>
#include <zlib.h>
#include <iostream>

class StreambufGZ : public std::streambuf {
private:
	static const int bufferSize = 47+256;    // size of data buff
	// totals 512 bytes under g++ for igzstream at the end.

	gzFile           file;               // file handle for compressed file
	char             buffer[bufferSize]; // data buffer
	char             opened;             // open/close state of stream
	int              mode;               // I/O mode

	int flush_buffer();
public:
	StreambufGZ() : file(0), opened(0), mode(0) {
		setp( buffer, buffer + (bufferSize-1));
		setg( buffer + 4, buffer + 4, buffer + 4);
		// ASSERT: both input & output capabilities will not be used together
	}
	int isOpen() { return opened; }
	StreambufGZ* open( const char* name, int open_mode);
	StreambufGZ* close();
	~StreambufGZ() { close(); }

	virtual int     overflow( int c = EOF);
	virtual int     underflow();
	virtual int     sync();
};

class StreamGZ_Base : virtual public std::ios {
protected:
	StreambufGZ buf;
public:
	StreamGZ_Base() { init(&buf); }
	StreamGZ_Base( const char* name, int open_mode);
	~StreamGZ_Base();

	void open( const char* name, int open_mode);
	void close();
	StreambufGZ* rdbuf() { return &buf; }
};

// ----------------------------------------------------------------------------
// User classes. Use igzstream and ogzstream analogously to ifstream and
// ofstream respectively. They read and write files based on the gz* 
// function interface of the zlib. Files are compatible with gzip compression.
// ----------------------------------------------------------------------------

class InputStreamGZ : public StreamGZ_Base, public std::istream {
public:
	InputStreamGZ() : std::istream( &buf) {}
	InputStreamGZ( const char* name, int open_mode = std::ios::in)
	: StreamGZ_Base( name, open_mode), std::istream( &buf) {}
	StreambufGZ* rdbuf() { return StreamGZ_Base::rdbuf(); }
	void open( const char* name, int open_mode = std::ios::in) {
		StreamGZ_Base::open( name, open_mode);
	}
};

class OutputStreamGZ : public StreamGZ_Base, public std::ostream {
public:
	OutputStreamGZ() : std::ostream( &buf) {}
	OutputStreamGZ( const char* name, int mode = std::ios::out)
	: StreamGZ_Base( name, mode), std::ostream( &buf) {}

	StreambufGZ* rdbuf() { return StreamGZ_Base::rdbuf(); }
	void open( const char* name, int open_mode = std::ios::out) {
		StreamGZ_Base::open( name, open_mode);
	}
};

#endif // STREAMGZ_HPP_
