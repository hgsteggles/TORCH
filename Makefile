# MAKEFILE FOR simple C++ programming

CFLAGS = -O3 -g -pedantic -Wall -std=c++11
INCLUDE = -I/home/harry/c++ -I/usr/local/include/
#LIBS =  -lreadline -lncurses
#LIBS = -lhdf5
#LIBS = /home/harry/Documents/libraries/cavlib/cavlib.a -lgsl -lgslcblas -lfftw3
CXX = g++
SRCDIR = src
HDRDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp
HDREXT = hpp
OBJEXT = o
SRCS = main.cpp io.cpp hydro.cpp grid3d.cpp gridcell.cpp parameters.cpp rtmodule.cpp
HDRS = $(SRCS:.$(SRCEXT)=.$(HDREXT))
OBJ = $(SRCS:.$(SRCEXT)=.$(OBJEXT))

%.o : $(SRCDIR)/%.$(SRCEXT) $(HDRDIR)/%.$(HDREXT)
	$(CXX) -c $(CFLAGS) $(INCLUDE) $< -o $(OBJDIR)/$@

main : $(OBJ)
	$(CXX) $(CFLAGS) $(INCLUDE) -o $(BINDIR)/$@ $(addprefix $(OBJDIR)/, $^)

.PHONY: clean
clean:
	rm obj/*.o bin/main
