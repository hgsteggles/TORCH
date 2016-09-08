# MAKEFILE FOR simple C++ programming

CFLAGS = -O2 -g -pedantic -Wall -std=c++11
INCLUDE = -I./src -I./src/Torch -I./src/MPI -I./src/IO -I./src/Fluid -I./src/Integrators -I./src/Misc -I./include/ -I./include/lua-5.2.3/
LIBS = -L./lib/ -l:liblua.a -lz -ldl
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
IO/DataReader \
IO/Logger \
IO/ProgressBar \
IO/StreamGZ \
Torch/Constants \
Torch/Converter \
Torch/Parameters \
Fluid/Fluid \
Fluid/GridCellCollection \
Fluid/Grid \
Fluid/GridCell \
Fluid/Star \
Fluid/PartitionManager \
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
	mkdir -p $(dir $@)
	$(MPICXX) -c $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS)
	
$(OBJDIR)/%.o : $(SRCDIR)/%.$(SRCEXT)
	mkdir -p $(dir $@)
	$(MPICXX) -c $(CFLAGS) $< -o $@ $(INCLUDE) $(LIBS)

torch : $(FULLPATHOBJ)
	$(MPICXX) $(CFLAGS) -o $@ $^ $(INCLUDE) $(LIBS)

test : $(FULLPATHOBJ)
	$(MPICXX) $(CFLAGS) -o $@ $^ $(INCLUDE) $(LIBS)

.PHONY: clean
clean:
	rm $(FULLPATHOBJ) torch test
