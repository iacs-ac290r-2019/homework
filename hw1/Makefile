# Include the configuration file
include config.mk

# Lists of intermediate object files (not executables)
objs=PoissonProblem.o GaussQuad.o

# Corresponding source files for these object files
src=$(patsubst %.o,%.cpp,$(objs) $(mg_objs))

# List of executables with .x suffix
execs_x= poisson.x

# Lists of executables (prefix only)
execs=$(patsubst %.x,%,$(execs_x) $(mg_execs))

# Define target for make all
all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs_x)

# Build object files from source
%.o: %.cpp
	$(cxx) $(includes) $(cflags) $(linkage) -c $<

# Build and Link executables from object files
poisson: poisson.cpp $(objs)
	$(cxx) $(includes) $(cflags) $(linkage) $(cxx_out)

# The PHONY command tells GNU Make than "clean" is a phony target, which means
# that it doesn't correspond to an actual file called "clean'.
.PHONY: clean depend
