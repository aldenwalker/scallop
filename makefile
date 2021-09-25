CC=g++
CFLAGS=-O3 -fcommon #-g -Wall
IFLAGS=-I/sw/include -I/opt/local/include
LDFLAGS=-L/sw/lib -I/opt/local/lib -lglpk -lgmp

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
	$(CC) $(CFLAGS) $(IFLAGS) -c rational.cc

word.o: word.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c word.cc

lp.o: lp.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c lp.cc

lp.o_GUR: lp.cc
	$(CC) $(CFLAGS) $(IFLAGS) $(GURINC) -DGUROBI_INSTALLED -c lp.cc

scallop.o: scallop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scallop.cc

scallop_with_gurobi: $(DIRS) scallop.o rational.o word.o lp.o_GUR scylla gallop trollop exlp-package
	$(CC) $(CFLAGS) -o scallop *.o exlp-package/*.o scylla/*.o gallop/*.o trollop/*.o scabble/*.o hallop/*.o $(GURLIB) $(LDFLAGS)

scallop: $(DIRS) scallop.o rational.o word.o lp.o scylla gallop trollop scabble exlp-package
	$(CC) $(CFLAGS) -o scallop *.o exlp-package/*.o scylla/*.o gallop/*.o trollop/*.o scabble/*.o hallop/*.o $(LDFLAGS)

clean: 
	rm *.o
	rm exlp-package/*.o
	rm scylla/*.o
	rm gallop/*.o
	rm trollop/*.o
	rm scabble/*.o
	rm hallop/*.o
	rm scallop
