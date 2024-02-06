TARGET=
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
	$(CC) $(CFLAGS) $(CXXFLAGS)  $(IFLAGS) -c rational.cc $(TARGET)

word.o: word.cc
	$(CC) $(CFLAGS) $(CXXFLAGS) $(IFLAGS) -c word.cc $(TARGET)

lp.o: lp.cc
	$(CC) $(CFLAGS) $(CXXFLAGS) $(IFLAGS) -c lp.cc $(TARGET)

lp.o_GUR: lp.cc
	$(CC) $(CFLAGS) $(CXXFLAGS) $(IFLAGS) $(GURINC) -DGUROBI_INSTALLED -c lp.cc $(TARGET)

scallop.o: scallop.cc
	$(CC) $(CFLAGS) $(CXXFLAGS) $(IFLAGS) -c scallop.cc $(TARGET)

scallop_with_gurobi: $(DIRS) scallop.o rational.o word.o lp.o_GUR scylla gallop trollop exlp-package
	$(CC) $(CFLAGS) $(CXXFLAGS) -o scallop *.o exlp-package/*.o scylla/*.o gallop/*.o trollop/*.o scabble/*.o hallop/*.o $(GURLIB) $(LDFLAGS) $(TARGET)

scallop: $(DIRS) scallop.o rational.o word.o lp.o scylla gallop trollop scabble exlp-package
	$(CC) $(CFLAGS) $(CXXFLAGS) -o scallop *.o exlp-package/*.o scylla/*.o gallop/*.o trollop/*.o scabble/*.o hallop/*.o $(LDFLAGS) $(TARGET)

clean:
	rm *.o
	rm exlp-package/*.o
	rm scylla/*.o
	rm gallop/*.o
	rm trollop/*.o
	rm scabble/*.o
	rm hallop/*.o
	rm scallop
