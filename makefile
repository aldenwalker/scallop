CC=g++
CFLAGS=-O3 -fcommon  -g -Wall 
CXXFLAGS=-std=c++11
IFLAGS=-I/usr/local/include -I${CONDA_PREFIX}/include
LDFLAGS=-I/usr/local/include  -L/usr/local/lib -lglpk -lgmpxx -lgmp 
 

#gurobi stuff
GURDIR = /home/akwalker/Documents/software/gurobi501/linux64
GURINCDIR = $(GURDIR)/include
GURLIBDIR = $(GURDIR)/lib
GURINC = -I$(GURINCDIR) 
GURLIB = -L$(GURLIBDIR) -lgurobi50

DIRS = exlp-package scylla gallop trollop scabble hallop

all: scallop

.PHONY : $(DIRS)
$(DIRS) :
	$(MAKE) -C $@

rational.o: rational.cc
	$(CC) $(CFLAGS) $(CXXFLAGS)  $(IFLAGS) -c rational.cc -target x86_64-apple-macos14

word.o: word.cc
	$(CC) $(CFLAGS) $(CXXFLAGS) $(IFLAGS) -c word.cc -target x86_64-apple-macos14

lp.o: lp.cc
	$(CC) $(CFLAGS) $(CXXFLAGS) $(IFLAGS) -c lp.cc -target x86_64-apple-macos14

lp.o_GUR: lp.cc
	$(CC) $(CFLAGS) $(CXXFLAGS) $(IFLAGS) $(GURINC) -DGUROBI_INSTALLED -c lp.cc -target x86_64-apple-macos14

scallop.o: scallop.cc
	$(CC) $(CFLAGS) $(CXXFLAGS) $(IFLAGS) -c scallop.cc -target x86_64-apple-macos14

scallop_with_gurobi: $(DIRS) scallop.o rational.o word.o lp.o_GUR scylla gallop trollop exlp-package
	$(CC) $(CFLAGS) $(CXXFLAGS) -o scallop *.o exlp-package/*.o scylla/*.o gallop/*.o trollop/*.o scabble/*.o hallop/*.o $(GURLIB) $(LDFLAGS) -target x86_64-apple-macos14

scallop: $(DIRS) scallop.o rational.o word.o lp.o scylla gallop trollop scabble exlp-package
	$(CC) $(CFLAGS) $(CXXFLAGS) -o scallop *.o exlp-package/*.o scylla/*.o gallop/*.o trollop/*.o scabble/*.o hallop/*.o $(LDFLAGS) -target x86_64-apple-macos14

clean: 
	rm *.o
	rm exlp-package/*.o
	rm scylla/*.o
	rm gallop/*.o
	rm trollop/*.o
	rm scabble/*.o
	rm hallop/*.o
	rm scallop
