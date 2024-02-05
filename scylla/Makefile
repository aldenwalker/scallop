CC=g++
CFLAGS=-O3 -std=c++14 #-g -Wall
IFLAGS=-I${CONDA_PREFIX}/include

all: scylla.o scylla_classes.o 

scylla.o: scylla.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scylla.cc -target x86_64-apple-macos14

scylla_classes.o: scylla_classes.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scylla_classes.cc -target x86_64-apple-macos14

clean: 
	rm *.o

