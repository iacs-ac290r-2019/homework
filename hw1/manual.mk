PoissonProblem.o : PoissonProblem.cpp
	g++ -c -Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include -o PoissonProblem.o PoissonProblem.cpp

poisson.o : poisson.cpp
	g++ -c -Wall -ansi -std=c++17 -O3 -I/d/SoftwareLibraries/yaml-cpp/include -o poisson.o poisson.cpp

poisson.x : poisson.cpp
	g++ -Wall -ansi -std=c++17 -O3 -o poisson.x poisson.o PoissonProblem.o libyaml-cpp.a

all: PoissonProblem.o poisson.o poisson.x
