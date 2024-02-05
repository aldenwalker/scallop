CC=gcc 
CFLAGS=-DHAVE_CONFIG_H -fcommon -I.  -I/usr/local/include -I${CONDA_PREFIX}/include -DNO_GMP_HASH -O3

all: main.o hash.o solve_ip.o preprocess.o matrix.o lpstruct.o \
     eta_file.o solve_lp.o mps.o lu.o EXLPvector.o mylib.o

EXLPvector.o : *.c  *.h
	$(CC) $(CFLAGS) -c EXLPvector.c  -target x86_64-apple-macos14
	
lu.o : *.c  *.h
	$(CC) $(CFLAGS) -c lu.c  -target x86_64-apple-macos14
	
mps.o : *.c  *.h
	$(CC) $(CFLAGS) -c mps.c  -target x86_64-apple-macos14
	
solve_lp.o : *.c  *.h
	$(CC) $(CFLAGS) -c solve_lp.c  -target x86_64-apple-macos14
	
eta_file.o : *.c  *.h
	$(CC) $(CFLAGS) -c eta_file.c  -target x86_64-apple-macos14
	
lpstruct.o : *.c  *.h
	$(CC) $(CFLAGS) -c lpstruct.c  -target x86_64-apple-macos14
	
matrix.o : *.c  *.h
	$(CC) $(CFLAGS) -c matrix.c  	-target x86_64-apple-macos14
	
preprocess.o : *.c  *.h
	$(CC) $(CFLAGS) -c preprocess.c  	-target x86_64-apple-macos14

solve_ip.o : *.c  *.h
	$(CC) $(CFLAGS) -c solve_ip.c  -target x86_64-apple-macos14

hash.o : *.c  *.h
	$(CC) $(CFLAGS) -c hash.c  	-target x86_64-apple-macos14

main.o : *.c  *.h
	$(CC) $(CFLAGS) -c main.c  -target x86_64-apple-macos14
	
mylib.o : *.c  *.h
	$(CC) $(CFLAGS) -c mylib.c -target x86_64-apple-macos14
	
clean:
	rm *.o
