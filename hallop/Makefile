CC=g++
CFLAGS=-g -Wall
IFLAGS=-I${CONDA_PREFIX}/include

all: hallop.o free_group_chain.o pieces.o hallop_lp.o

hallop.o: hallop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c hallop.cc -target x86_64-apple-macos14

free_group_chain.o: free_group_chain.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c free_group_chain.cc -target x86_64-apple-macos14

pieces.o: pieces.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c pieces.cc -target x86_64-apple-macos14

hallop_lp.o: hallop_lp.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c hallop_lp.cc -target x86_64-apple-macos14
clean: 
	rm *.o

