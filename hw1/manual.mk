# Compile PoissonProblem class
PoissonProblem.o : PoissonProblem.cpp
	g++ -c -Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include -o PoissonProblem.o PoissonProblem.cpp

# Compile poisson program
poisson.o : poisson.cpp
	g++ -c -Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include -o poisson.o poisson.cpp

# Link poisson program
poisson.x : poisson.o PoissonProblem.o
	g++ -Wall -ansi -std=c++17 -O3 -o poisson.x poisson.o PoissonProblem.o -L /d/SoftwareLibraries/lib -lyaml-cpp

# Targets
all: PoissonProblem.o poisson.o poisson.x
poisson: PoissonProblem.o poisson.o poisson.x

# A Makefile target to remove all the built files
clean:
	rm -f PoissonProblem.o poisson.o poisson.x

.PHONY: clean depend
