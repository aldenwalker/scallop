CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/sw/include

all: hallop.o free_group_chain.o 

hallop.o: hallop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c hallop.cc

free_group_chain.o: free_group_chain.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c free_group_chain.cc

clean: 
	rm *.o

