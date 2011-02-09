CC=g++
CFLAGS=-g -Wall
IFLAGS=-I/sw/include
LDFLAGS=-L/sw/lib -lglpk -lgmp
all: scallop

scallop.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c scallop.cc
	
word.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c word.cc
	
draw.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c draw.cc

rational.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c rational.cc
	
lp.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c lp.cc
	
io.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c io.cc


.PHONY: exlp-package
exlp-package: 
	cd exlp-package; make

scallop: exlp-package scallop.o word.o draw.o rational.o lp.o io.o
	$(CC) $(CFLAGS) -o scallop scallop.o word.o draw.o rational.o lp.o io.o exlp-package/*.o $(LDFLAGS)
	
	
clean: 
	rm scallop
	rm *.o
	cd exlp-package; rm *.o
