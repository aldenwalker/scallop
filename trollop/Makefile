CC=g++
CFLAGS=-O3 #-g -Wall
IFLAGS=-I${CONDA_PREFIX}/include

all: trollop.o trollop_classes.o

trollop.o: trollop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c trollop.cc -target x86_64-apple-macos14

trollop_classes.o: trollop_classes.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c trollop_classes.cc -target x86_64-apple-macos14

clean: 
	rm *.o


