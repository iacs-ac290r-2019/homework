# C++ compiler
# cxx=g++ -fopenmp
cxx=g++

# Compilation flags
cflags=-Wall -ansi -std=c++17 -O3

# Output command
# Note $@ is a shorthand for the file to be built and $^ is a shorthand for all the dependencies. 
# (Can also use $< to mean the first dependency only.)
cxx_out=-o $@.x $^

# BLAS/LAPACK flags for linear algebra
lapack_link=-llapack -lblas

# Root directory for manually installed software libraries; environment variable SOFTWARE_LIBRARY_DIR
software_libs_include=$(addprefix -I,$(SOFTWARE_LIBRARY_DIR))
software_libs_link=$(addprefix -L, $(addsuffix /lib, $(SOFTWARE_LIBRARY_DIR)))

# Boost flags are read using environment variable BOOST_DIR
boost_include=$(addprefix -I,$(BOOST_DIR))

# yaml-cpp flags for parsing YAML configuration file
yaml_include=$(addsuffix /yaml-cpp/include, $(software_libs_include))
yaml_link=-lyaml-cpp

# Additional include directories
includes=$(boost_include) $(yaml_include)

# Linker flags
linkage=$(lapack_link) $(software_libs_link) $(yaml_link)
