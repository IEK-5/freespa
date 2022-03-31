SRC=freespa.c testspa.c spa.c
OBJ=freespa.o spa.o
HDR=freespa.h freespa_tables.h spa.h
CC=gcc
CFLAGS=-O3 -flto
LFLAGS=-lm

test: $(SRC) $(HDR) $(OBJ)
	$(CC) $(CFLAGS) -o testspa testspa.c $(OBJ) $(LFLAGS)
clean:
	rm *.o testspa
