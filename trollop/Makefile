CC=g++
CFLAGS=-O3 #-g -Wall
IFLAGS=-I/sw/include -I${CONDA_PREFIX}/include

all: trollop.o trollop_classes.o

trollop.o: trollop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c trollop.cc

trollop_classes.o: trollop_classes.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c trollop_classes.cc

clean: 
	rm *.o


