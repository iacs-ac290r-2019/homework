# *****************************************************************************
# Set environment
# Warn on unset variables
MAKEFLAGS=--warn-undefined-variables -j32

# Include the configuration file
include config.mk

# Include the file dependencies
include Makefile.dep

# *****************************************************************************
# Source, object and executable files
# List of source files for non-executables (e.g. classes); filename.cpp
SRC=PoissonProblem.cpp GaussQuad.cpp

# Lists of intermediate object files (not executables); filename.o
OBJ := $(patsubst %.cpp, %.o, $(SRC))

# Lists of executables (prefix only)
EXEC := poisson

# List of executables (with .x suffix)
EXEC_X=$(patsubst %, %.x, $(EXEC))

# *****************************************************************************
# Define target for make all
all: $(EXEC)

# A Makefile target to remove all the built files
clean:
	rm -f $(OBJ) $(EXEC_X)

# *****************************************************************************
# Build object files from source (not executables)
%.o: %.cpp
	$(cxx) -c $< \
		$(cflags) \
		$(includes) \


# *****************************************************************************
# Build and link executables from source and object files
poisson: poisson.cpp $(OBJ)
	$(cxx) $(cxx_out) \
		$(cflags) \
		$(includes) \
		$(linkage) \


# *****************************************************************************
# Set default goal to all
.DEFAULT_GOAL := all

# The PHONY command tells GNU Make than "clean" is a phony target, which means
# that it doesn't correspond to an actual file called "clean'.
.PHONY: clean
