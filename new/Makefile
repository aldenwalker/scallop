CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/sw/include
LDFLAGS=-L/sw/lib -lglpk -lgmp

#gurobi stuff
GURDIR = /home/awalker/Documents/research/software/gurobi500/linux64
GURINCDIR = $(GURDIR)/include
GURLIBDIR = $(GURDIR)/lib
GURINC = -I$(GURINCDIR) 
GURLIB = -L$(GURLIBDIR) -lgurobi50

all: scallop_with_gurobi

.PHONY: exlp-package
exlp-package: 
	cd exlp-package; make

rational.o: rational.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c rational.cc

lp.o: lp.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c LP.cc

scallop_with_gurobi: exlp-package rational.o lp.o


scallop: exlp-package scallop.o word.o draw.o rational.o lp.o io.o
	$(CC) $(CFLAGS) -o scallop scallop.o word.o draw.o rational.o lp.o io.o exlp-package/*.o $(LDFLAGS)

scylla: exlp-package scylla.o rational.o scylla_lp.o scylla_classes.o word.o
	$(CC) $(CFLAGS) -o scylla scylla.o scylla_lp.o scylla_classes.o rational.o word.o exlp-package/*.o $(CLNFLAGS) $(GURLIB) $(LDFLAGS)


clean: 
	rm *.o
	cd exlp-package; rm *.o
