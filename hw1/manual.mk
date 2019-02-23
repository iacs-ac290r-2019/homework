# Load dependencies
include Makefile.dep

# Compile PoissonFiniteElement class
PoissonFiniteElement.o : PoissonFiniteElement.hpp PoissonFiniteElement.cpp
	g++ -c PoissonFiniteElement.cpp -o PoissonFiniteElement.o \
	-Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include 

# Compile LinearSolve
LinearSolve.o : LinearSolve.hpp LinearSolve.cpp
	g++ -c LinearSolve.cpp -o LinearSolve.o \
	-Wall -ansi -std=c++17 -O3

# Compile PoissonExamples
PoissonExamples.o : PoissonExamples.hpp PoissonExamples.cpp
	g++ -c PoissonExamples.cpp -o PoissonExamples.o \
	-Wall -ansi -std=c++17 -O3

# Compile poisson program
poisson.o : poisson.cpp
	g++ -c poisson.cpp -o poisson.o \
	-Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include 

# Link poisson program
poisson.x : poisson.o PoissonFiniteElement.o PoissonExamples.o LinearSolve.o
	g++ -o poisson.x poisson.o PoissonFiniteElement.o PoissonExamples.o LinearSolve.o \
	-L /d/SoftwareLibraries/lib -llapack -lblas -lyaml-cpp

# Targets
all: poisson.x
poisson: poisson.x

# Set default goal to all
.DEFAULT_GOAL := all

# A Makefile target to remove all the built files
clean:
	rm -f PoissonFiniteElement.o PoissonExamples.o LinearSolve.o poisson.o poisson.x

.PHONY: clean depend
