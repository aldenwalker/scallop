CC=g++
CFLAGS=-g -Wall
IFLAGS=-I${CONDA_PREFIX}/include

all: scabble.o 

scabble.o: scabble.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scabble.cc -target x86_64-apple-macos14

clean: 
	rm *.o

