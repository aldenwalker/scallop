CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/sw/include

all: hallop.o free_group_chain.o pieces.o hallop_lp.o

hallop.o: hallop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c hallop.cc

free_group_chain.o: free_group_chain.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c free_group_chain.cc

pieces.o: pieces.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c pieces.cc

hallop_lp.o: hallop_lp.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c hallop_lp.cc
clean: 
	rm *.o

