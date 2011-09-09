CC=g++
CFLAGS=-g -pg -Wall
IFLAGS=-I/sw/include
LDFLAGS=-L/sw/lib -lglpk -lgmp
all: scallop scyllop

scallop.o: scallop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scallop.cc

word.o: word.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c word.cc

draw.o: draw.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c draw.cc

rational.o: rational.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c rational.cc

lp.o: lp.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c lp.cc

io.o: io.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c io.cc

scyllop.o: scyllop.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scyllop.cc
	
scyllop_lp.o:	scyllop_lp.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scyllop_lp.cc
	
scyllop_classes.o: scyllop_classes.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scyllop_classes.cc

.PHONY: exlp-package
exlp-package: 
	cd exlp-package; make

scallop: exlp-package scallop.o word.o draw.o rational.o lp.o io.o
	$(CC) $(CFLAGS) -o scallop scallop.o word.o draw.o rational.o lp.o io.o exlp-package/*.o $(LDFLAGS)

scyllop: exlp-package scyllop.o rational.o scyllop_lp.o scyllop_classes.o
	$(CC) $(CFLAGS) -o scyllop scyllop.o scyllop_lp.o scyllop_classes.o rational.o exlp-package/*.o $(LDFLAGS)


clean: 
	rm scallop
	rm *.o
	cd exlp-package; rm *.o
