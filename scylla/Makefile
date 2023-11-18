CC=g++
CFLAGS=-O3 #-g -Wall
IFLAGS=-I/sw/include -I${CONDA_PREFIX}/include

all: scylla.o scylla_classes.o 

scylla.o: scylla.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scylla.cc

scylla_classes.o: scylla_classes.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scylla_classes.cc

clean: 
	rm *.o

