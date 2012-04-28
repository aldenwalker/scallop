CC=g++
CFLAGS=-O3 #-g -Wall
IFLAGS=-I/sw/include
LDFLAGS=-L/sw/lib -lglpk -lgmp

#cplex stuff
SYSTEM     = x86-64_sles10_4.1
LIBFORMAT  = static_pic
CPLEXDIR=/opt/ibm/ILOG/CPLEX_Studio124/cplex
CPLEXINCDIR=$(CPLEXDIR)/include
CPLEXLIBDIR=$(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CLNFLAGS=-L$(CPLEXLIBDIR) -lcplex -lm -pthread
CPLEXFLAGS=-I$(CPLEXINCDIR)

all: scallop scylla

scallop.o: scallop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scallop.cc

word.o: word.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c word.cc

draw.o: draw.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c draw.cc

rational.o: rational.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c rational.cc

lp.o: lp.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c lp.cc

io.o: io.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c io.cc

scylla.o: scylla.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scylla.cc

scylla_lp.o:	scylla_lp.cc
	$(CC) $(CFLAGS) $(IFLAGS) $(CPLEXFLAGS) -c scylla_lp.cc

scylla_classes.o: scylla_classes.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scylla_classes.cc

.PHONY: exlp-package
exlp-package: 
	cd exlp-package; make

scallop: exlp-package scallop.o word.o draw.o rational.o lp.o io.o
	$(CC) $(CFLAGS) -o scallop scallop.o word.o draw.o rational.o lp.o io.o exlp-package/*.o $(LDFLAGS)

scylla: exlp-package scylla.o rational.o scylla_lp.o scylla_classes.o word.o
	$(CC) $(CFLAGS) -o scylla scylla.o scylla_lp.o scylla_classes.o rational.o word.o exlp-package/*.o $(CLNFLAGS) $(LDFLAGS)

clean: 
	rm scylla
	rm scallop
	rm *.o
	cd exlp-package; rm *.o
