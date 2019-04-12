# Detect the operating system
ifeq ($(OS),Windows_NT) 
    detected_OS := Windows
else
    detected_OS := $(shell sh -c 'uname 2>/dev/null || echo Unknown')
endif

# Show the operating system
# $(info OS is $(detected_OS))

# Set the compiler and flags based on the operating system
# Supported OS are Mac OSX ("Darwin") and Linux ("Linux")

# Is the OS Mac?
ifeq ($(detected_OS), Darwin)
cxx := g++-mp-7
cflags := -fopenmp -Wall -ansi -pedantic -O3
lp_lflags := -llapack -lblas
endif

# Is the OS Linux?
ifeq ($(detected_OS), Linux)
cxx := g++-8
cflags := -fopenmp -Wall -ansi -pedantic -O3
lp_lflags := -llapack -lblas
endif

# Lists of files to be built
objs=lbm.o file.o
src=$(patsubst %.o,%.cc,$(objs))
execs=lbm_pbc lbm_mbc

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

lbm_pbc: lbm_pbc.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

lbm_mbc: lbm_mbc.cc $(objs)
	$(cxx) $(cflags) -o $@ $^

.PHONY: clean depend
