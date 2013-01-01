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

DIRS = scylla exlp-package

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

scallop_with_gurobi: $(DIRS) scallop.o rational.o word.o lp.o_GUR exlp-package
	$(CC) $(CFLAGS) -o scallop *.o exlp-package/*.o scylla/*.o $(GURLIB) $(LDFLAGS)

scallop: $(DIRS) scallop.o rational.o word.o lp.o exlp-package
	$(CC) $(CFLAGS) -o scallop *.o exlp-package/*.o scylla/*.o $(LDFLAGS)

clean: 
	rm *.o
	cd exlp-package; rm *.o
