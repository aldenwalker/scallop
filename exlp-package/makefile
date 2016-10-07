CC=gcc 
CFLAGS=-DHAVE_CONFIG_H -I. -I/sw/include  -I/opt/local/include -DNO_GMP_HASH -O3

all: main.o hash.o solve_ip.o preprocess.o matrix.o lpstruct.o \
     eta_file.o solve_lp.o mps.o lu.o EXLPvector.o mylib.o

EXLPvector.o : *.c  *.h
	$(CC) $(CFLAGS) -c EXLPvector.c  
	
lu.o : *.c  *.h
	$(CC) $(CFLAGS) -c lu.c  
	
mps.o : *.c  *.h
	$(CC) $(CFLAGS) -c mps.c  
	
solve_lp.o : *.c  *.h
	$(CC) $(CFLAGS) -c solve_lp.c  
	
eta_file.o : *.c  *.h
	$(CC) $(CFLAGS) -c eta_file.c  
	
lpstruct.o : *.c  *.h
	$(CC) $(CFLAGS) -c lpstruct.c  
	
matrix.o : *.c  *.h
	$(CC) $(CFLAGS) -c matrix.c  	
	
preprocess.o : *.c  *.h
	$(CC) $(CFLAGS) -c preprocess.c  	

solve_ip.o : *.c  *.h
	$(CC) $(CFLAGS) -c solve_ip.c  

hash.o : *.c  *.h
	$(CC) $(CFLAGS) -c hash.c  	

main.o : *.c  *.h
	$(CC) $(CFLAGS) -c main.c  
	
mylib.o : *.c  *.h
	$(CC) $(CFLAGS) -c mylib.c 
	
clean:
	rm *.o
