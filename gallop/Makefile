CC=g++
CFLAGS=-O3 #-g -Wall
IFLAGS=-I/sw/include -I${CONDA_PREFIX}/include

all: gallop.o graph.o

gallop.o: gallop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c gallop.cc

graph.o: graph.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c graph.cc

clean: 
	rm *.o


