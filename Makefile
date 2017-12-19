shell = /bin/sh
CC = gcc -std=gnu99 -m64 -O3 -funroll-loops 

LIBS = -L/opt/local/lib 
FLAGS = -lgsl -lm

default: program

.SUFFIXES: .c .o
.c:
	$(CC) -c -o $*.o $<
.o:
	$(CC) -c -o $*.o $<

source:
	$(CC) ${LIBS} -c quickhull.c

program: quickhull.o
	$(CC) ${LIBS} ${FLAGS} -o hull.x *.o

clean:
	rm -f *.o *.x a.out
