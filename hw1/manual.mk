# Load dependencies
include Makefile.dep

# Compile PoissonProblem class
PoissonProblem.o : PoissonProblem.hpp PoissonProblem.cpp
	g++ -c PoissonProblem.cpp -o PoissonProblem.o \
		-Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include 

# Compile poisson program
poisson.o : poisson.cpp
	g++ -c poisson.cpp -o poisson.o \
		-Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include 

# Link poisson program
poisson.x : poisson.o PoissonProblem.o
	g++ -o poisson.x poisson.o PoissonProblem.o \
		-L /d/SoftwareLibraries/lib -lyaml-cpp

# Targets
all: poisson.x
poisson: poisson.x

# Set default goal to all
.DEFAULT_GOAL := all

# A Makefile target to remove all the built files
clean:
	rm -f PoissonProblem.o poisson.o poisson.x

.PHONY: clean depend
