# Load dependencies
include Makefile.dep

# Compile PoissonFiniteElement class
PoissonFiniteElement.o : PoissonFiniteElement.hpp PoissonFiniteElement.cpp
	g++ -c PoissonFiniteElement.cpp -o PoissonFiniteElement.o \
		-Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include 

# Compile poisson program
poisson.o : poisson.cpp
	g++ -c poisson.cpp -o poisson.o \
		-Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include 

# Link poisson program
poisson.x : poisson.o PoissonFiniteElement.o
	g++ -o poisson.x poisson.o PoissonFiniteElement.o \
		-L /d/SoftwareLibraries/lib -lyaml-cpp

# Targets
all: poisson.x
poisson: poisson.x

# Set default goal to all
.DEFAULT_GOAL := all

# A Makefile target to remove all the built files
clean:
	rm -f PoissonFiniteElement.o poisson.o poisson.x

.PHONY: clean depend
