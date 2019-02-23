# *************************************************************************************************
# Newline and tab characters for legibile printing to screen
# *************************************************************************************************
# newline character
define NEWLINE


endef

# tab character
TAB := $(shell printf '\011')

# *************************************************************************************************
# Make settings, environment variables, initial status message
# *************************************************************************************************
# Get the number of jobs for make from the environment; default to all of them
ifndef MAKE_JOBS
	NUMPROC := $(shell grep -c ^processor /proc/cpuinfo)
	MAKE_JOBS := $(shell echo $(NUMPROC)/1 + 0 | bc)
endif

# Make settings: warn on unset variables, use parallel processing
MAKEFLAGS=--warn-undefined-variables --jobs=$(MAKE_JOBS)

# Show the number of jobs
$(info Running make with up to $(MAKE_JOBS) parallel jobs.$(NEWLINE))

# *************************************************************************************************
# Compiler settings
# *************************************************************************************************
# C++ compiler
CXX=g++

# Compilation flags
CXX_FLAGS= \$(NEWLINE) $(TAB) -fopenmp -Wall -ansi -std=c++17 -O3

# Output command for object files
# Note $< is a shorthand for the first dependency
# CXX_OUT_OBJ = -c $<
CXX_OUT_OBJ = -c $< -o $@

# Output command for executables
# Note $@ is a shorthand for the file to be built
# Note $^ is a shorthand for all the dependencies
CXX_OUT_EXE = -o $@.x $^

# *************************************************************************************************
# INCLUDE: Include paths for additional software libraries 
# *************************************************************************************************
# Root directory for manually installed software libraries; 
# found in environment variable SOFTWARE_LIBRARY_DIR
INCLUDE_USR := $(addprefix -I,$(SOFTWARE_LIBRARY_DIR))

# Boost flags are read using environment variable BOOST_DIR
ifdef BOOST_DIR
	INCLUDE_BOOST=$(addprefix -I,$(BOOST_DIR))
endif

# yaml-cpp flags for parsing YAML configuration file
INCLUDE_YAML=$(addsuffix /yaml-cpp/include, $(INCLUDE_USR))

# *************************************************************************************************
# LD: Linker flags for additional libraries
# *************************************************************************************************
# Directory with library files (.a) for additional software libraries
LD_FLAGS_USR := $(addprefix -L, $(addsuffix /lib, $(SOFTWARE_LIBRARY_DIR)))

# Math library
LD_LIB_MATH := -lm

# BLAS/LAPACK for linear algebra
LD_LIB_LAPACK := -llapack -lblas

# yaml-cpp for YAML file parsing
LD_LIB_YAML := -lyaml-cpp

# *************************************************************************************************
# Generate command line arguments 
# -I: additional include directories
# -L: additional library directories
# -l: additional libraries
# *************************************************************************************************
# Additional include directories -I
# Initialize INCLUDE to an empty string
INCLUDE := 
# Only include boost if it was supplied manually
ifdef BOOST_DIR
	INCLUDE := \$(NEWLINE) $(TAB) $(INCLUDE_BOOST)
endif

# yaml-cpp
INCLUDE := $(INCLUDE)\
\$(NEWLINE) $(TAB) $(INCLUDE_YAML)

# Additional library directories -L
LD_FLAGS := \
\$(NEWLINE) $(TAB) $(LD_FLAGS_USR) 

# Additional libraries -l
LD_LIBS := \
\$(NEWLINE) $(TAB) $(LD_LIB_MATH) $(LD_LIB_LAPACK) $(LD_LIB_YAML)
