CC=g++
CFLAGS=-O3 #-g -Wall
IFLAGS=-I${CONDA_PREFIX}/include

all: gallop.o graph.o

gallop.o: gallop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c gallop.cc -target x86_64-apple-macos14

graph.o: graph.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c graph.cc -target x86_64-apple-macos14

clean: 
	rm *.o


