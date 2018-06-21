TARGET = quickhull
CC = gcc 
CFLAGS = -std=gnu99 -m64 -O3 -funroll-loops -Wall -I/opt/local/include
LIBS = -lm -lgsl -L/opt/local/lib 

.PHONY: default all clean

default: $(TARGET)

all: default

.SUFFIXES: .c .o
.o .c:
	$(CC) $(CFLAGS) -c $< -o $@

source:
	$(CC) ${LIBS} ${INCLUDE} -c quickhull.c main.c

quickhull: quickhull.o
	$(CC) ${LIBS} ${FLAGS} -o quick.x quickhull.o

giftwrap: main.o
	$(CC) ${LIBS} ${FLAGS} -o wrap.x main.o


clean:
	rm -f *.o *.x a.out
