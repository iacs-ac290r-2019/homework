# Warn on unset variables
MAKEFLAGS=--warn-undefined-variables -j32

# Include the configuration file
include config.mk

# List of source files for non-executables (e.g. classes); filename.cpp
src=PoissonProblem.cpp GaussQuad.cpp

# Lists of intermediate object files (not executables); filename.o
objs = $(patsubst %.cpp,%.o,$(src))

# List of executables with .x suffix
execs_x = poisson.x

# Lists of executables (prefix only)
execs=$(patsubst %.x,%,$(execs_x))

# Define target for make all
all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs_x)

# Build object files from source, including executables
%.o: %.cpp
	$(cxx) -c $< \
		$(cflags) \
		$(includes) 

# Build and link executables
poisson: poisson.cpp $(objs)
	$(cxx) $(cxx_out) \
		$(cflags) \
		$(includes) \
		$(linkage) 

# The PHONY command tells GNU Make than "clean" is a phony target, which means
# that it doesn't correspond to an actual file called "clean'.
.PHONY: clean depend
