CC=gcc
CFLAGS=-DHAVE_CONFIG_H -fcommon -I.  -I/usr/local/include -I${CONDA_PREFIX}/include -DNO_GMP_HASH -O3

all: main.o hash.o solve_ip.o preprocess.o matrix.o lpstruct.o \
     eta_file.o solve_lp.o mps.o lu.o EXLPvector.o mylib.o

EXLPvector.o : *.c  *.h
	$(CC) $(CFLAGS) -c EXLPvector.c  $(TARGET)

lu.o : *.c  *.h
	$(CC) $(CFLAGS) -c lu.c  $(TARGET)

mps.o : *.c  *.h
	$(CC) $(CFLAGS) -c mps.c  $(TARGET)

solve_lp.o : *.c  *.h
	$(CC) $(CFLAGS) -c solve_lp.c  $(TARGET)

eta_file.o : *.c  *.h
	$(CC) $(CFLAGS) -c eta_file.c  $(TARGET)

lpstruct.o : *.c  *.h
	$(CC) $(CFLAGS) -c lpstruct.c  $(TARGET)

matrix.o : *.c  *.h
	$(CC) $(CFLAGS) -c matrix.c  	$(TARGET)

preprocess.o : *.c  *.h
	$(CC) $(CFLAGS) -c preprocess.c  	$(TARGET)

solve_ip.o : *.c  *.h
	$(CC) $(CFLAGS) -c solve_ip.c  $(TARGET)

hash.o : *.c  *.h
	$(CC) $(CFLAGS) -c hash.c  	$(TARGET)

main.o : *.c  *.h
	$(CC) $(CFLAGS) -c main.c  $(TARGET)

mylib.o : *.c  *.h
	$(CC) $(CFLAGS) -c mylib.c $(TARGET)

clean:
	rm *.o
