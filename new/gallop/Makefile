CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/sw/include

all: gallop.o graph.o

gallop.o: gallop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c gallop.cc

graph.o: graph.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c graph.cc

clean: 
	rm *.o


