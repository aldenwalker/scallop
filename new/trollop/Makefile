CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/sw/include

all: trollop.o trollop_classes.o

trollop.o: trollop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c trollop.cc

trollop_classes.o: trollop_classes.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c trollop_classes.cc

clean: 
	rm *.o


