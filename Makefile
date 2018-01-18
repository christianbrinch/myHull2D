shell = /bin/sh
CC = gcc -std=gnu99 -m64 -O3 -funroll-loops 

LIBS = -L/opt/local/lib 
FLAGS = -lgsl -lm

default: quickhull

.SUFFIXES: .c .o
.c:
	$(CC) -c -o $*.o $<
.o:
	$(CC) -c -o $*.o $<

source:
	$(CC) ${LIBS} -c quickhull.c main.c

quickhull: quickhull.o
	$(CC) ${LIBS} ${FLAGS} -o quick.x quickhull.o

giftwrap: main.o
	$(CC) ${LIBS} ${FLAGS} -o wrap.x main.o


clean:
	rm -f *.o *.x a.out
