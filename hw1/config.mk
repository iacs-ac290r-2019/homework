# C++ compiler
cxx=g++ -fopenmp

# Compilation flags
cflags=-Wall -ansi -std=c++17 -O3

# Output command
# Note $@ is a shorthand for the file to be built and $^ is a shorthand for all the dependencies. 
# (Can also use $< to mean the first dependency only.)
cxx_out=-o $@.x $^

# BLAS/LAPACK flags for linear algebra
lp_link=-llapack -lblas

# Boost flags are read using environment variable BOOST_DIR
boost_include=$(addprefix -I,$(BOOST_DIR))
boost_link=$(addsuffix($(addprefix -l,$(BOOST_DIR)),"/bin.v2/libs")

# yaml-cpp flags for parsing YAML configuration file
# TODO

# Additional include directories
includes=$(boost_include)

# Linker flags
linkage=$(lp_link)
