CC=g++
CFLAGS=-O3 -std=c++14 #-g -Wall
IFLAGS=-I${CONDA_PREFIX}/include

all: scylla.o scylla_classes.o

scylla.o: scylla.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scylla.cc $(TARGET)

scylla_classes.o: scylla_classes.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scylla_classes.cc $(TARGET)

clean:
	rm *.o

