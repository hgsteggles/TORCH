# MAKEFILE FOR simple C++ programming

CFLAGS = -O2 -g -pedantic -Wall -std=c++11
INCLUDE = -I./src -I./src/Torch -I./src/MPI -I./src/IO -I./src/Fluid -I./src/Integrators -I./src/Misc -I./include/ -I./include/lua-5.2.3/
LIBS = -L./lib/ -l:liblua.a -lz
CXX = mpic++
MPICXX = mpic++
SRCDIR = src
HDRDIR = src
OBJDIR = obj
BINDIR = bin
SRCEXT = cpp
HDREXT = hpp
OBJEXT = o

SUBDIRS = Torch MPI IO Fluid Integrators Misc

$(shell mkdir -p $(OBJDIR))
$(shell mkdir -p $(addprefix $(OBJDIR)/, $(SUBDIRS)))

FILES = main \
Torch/Torch \
MPI/MPI_Wrapper \
IO/DataPrinter \
IO/Logger \
IO/ProgressBar \
IO/StreamGZ \
Torch/Constants \
Torch/Converter \
Torch/Parameters \
Fluid/Fluid \
Fluid/GridFactory \
Fluid/Grid \
Fluid/Boundary \
Fluid/GridCell \
Fluid/Star \
Integrators/Hydro \
Integrators/Riemann \
Integrators/SlopeLimiter \
Integrators/Radiation \
Integrators/Thermodynamics \
Integrators/SplineData \
Misc/Timer
				
SRCS = $(FILES:=.$(SRCEXT))
HDRS = $(FILES:=.$(HDREXT))
OBJS = $(FILES:=.$(OBJEXT))
#HDRS = $(SRCS:.$(SRCEXT)=.$(HDREXT))
#OBJS = $(SRCS:.$(SRCEXT)=.$(OBJEXT))
FULLPATHOBJ = $(addprefix $(OBJDIR)/, $(OBJS))
FULLPATHSRC = $(addprefix $(SRCDIR)/, $(SRCS))
FULLPATHHDR = $(addprefix $(HDRDIR)/, $(HDRS))

$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT) $(HDRDIR)/%.$(HDREXT)
	@echo -n "Linking:   "
	$(CXX) -c $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS)
	
$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT)
	@echo -n "Linking:   "
	$(CXX) -c $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS)

torch : $(FULLPATHOBJ)
	@echo -n "Compiling: "
	$(CXX) $(CFLAGS) -o $@ $^ $(INCLUDE) $(LIBS)

test : $(FULLPATHOBJ)
	@echo -n "Compiling: "
	$(CXX) $(CFLAGS) -o $@ $^ $(INCLUDE) $(LIBS)

.PHONY: clean
clean:
	rm $(FULLPATHOBJ) torch test
