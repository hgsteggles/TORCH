# MAKEFILE FOR simple C++ programming

CFLAGS = -O3 -g -pedantic -Wall -std=c++11
#INCLUDE = -I/home/harry/c++ -I/usr/local/include/
INCLUDE = -I./libs
#LIBS =  -lreadline -lncurses
#LIBS = -lhdf5
#LIBS = /home/harry/Documents/libraries/cavlib/cavlib.a -lgsl -lgslcblas -lfftw3
MPILIBS = -lboost_mpi -lboost_serialization -lboost_system -lboost_filesystem -lboost_graph_parallel -lboost_iostreams /usr/local/lib/libprofiler.so
CXX = mpic++
MPICXX = mpic++
SRCDIR = src
HDRDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp
HDREXT = hpp
OBJEXT = o
FILES = main integrator mpihandler io parameters hydro radiation thermodynamics grid3d boundary external partition gridcell slopelimiter recipes
SRCS = $(FILES:=.$(SRCEXT))
HDRS = $(FILES:=.$(HDREXT))
OBJS = $(FILES:=.$(OBJEXT))
#HDRS = $(SRCS:.$(SRCEXT)=.$(HDREXT))
#OBJS = $(SRCS:.$(SRCEXT)=.$(OBJEXT))
FULLPATHOBJ = $(addprefix $(OBJDIR)/, $(OBJS))
FULLPATHSRC = $(addprefix $(SRCDIR)/, $(SRCS))
FULLPATHHDR = $(addprefix $(HDRDIR)/, $(HDRS))

$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(HDRDIR)/%.$(HDREXT)
	@echo -n "Linking: "
	$(CXX) -c $(CFLAGS) $< -o $@ $(INCLUDE)
	
$(OBJDIR)/mpihandler2.o : $(SRCDIR)/mpihandler.cpp $(HDRDIR)/mpihandler.hpp
	$(CXX) -c $(CFLAGS) $< -o $@ $(INCLUDE) $(MPILIBS)

main : $(FULLPATHOBJ)
	@echo -n "Compiling: "
	$(CXX) $(CFLAGS) -o $@ $^ $(INCLUDE) $(MPILIBS)

.PHONY: clean
clean:
	rm obj/*.o main
