CC=g++
CFLAGS=-g -Wall
IFLAGS=-I${CONDA_PREFIX}/include

all: scabble.o

scabble.o: scabble.cc
	$(CC) $(CFLAGS) $(IFLAGS) -c scabble.cc $(TARGET)

clean:
	rm *.o

