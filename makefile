CC=g++
CFLAGS=-O3 #-g -Wall
LDFLAGS=-lglpk -lgmp
all: scallop

exlp:
	cd exlp-package; make

scallop: exlp scallop.cc word.cc draw.cc rational.cc lp.cc
	$(CC) $(CFLAGS) -c scallop.cc word.cc draw.cc rational.cc lp.cc 
	$(CC) $(CFLGAS) -o scallop scallop.o word.o draw.o rational.o lp.o exlp-package/*.o $(LDFLAGS)
	
	
clean: 
	rm scallop
	rm *.o
	cd exlp-package; rm *.o
