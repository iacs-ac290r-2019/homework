# Set the compiler and flags
cxx=g++-mp-7 -O3 
cflags=-fopenmp -Wall -ansi -pedantic
lp_lflags=-llapack -lblas

# Lists of files to be built
objs=lbm.o file.o
src=$(patsubst %.o,%.cc,$(objs))
execs=lbm_test

all: $(execs)

# Include the file dependencies
-include Makefile.dep

# A Makefile target to refresh the dependency file
depend:
	$(cxx) -MM $(src) >Makefile.dep

# A Makefile target to remove all the built files
clean:
	rm -f $(objs) $(execs)

%.o: %.cc
	$(cxx) $(cflags) -c $<

lbm_test: lbm_test.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend
