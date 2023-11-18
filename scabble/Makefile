CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/sw/include -I${CONDA_PREFIX}/include

all: scabble.o 

scabble.o: scabble.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scabble.cc

clean: 
	rm *.o

