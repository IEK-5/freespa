SRC=freespa.c testspa.c
OBJ=freespa.o
HDR=freespa.h freespa_tables.h freespa_dt_table.h
CC=gcc
CFLAGS=-O3 -flto
LFLAGS=-lm


test: $(SRC) $(HDR) $(OBJ)
	$(CC) $(CFLAGS) -o testspa testspa.c $(OBJ) $(LFLAGS)

compare: $(SRC) $(HDR) $(OBJ) spa.c spa.h spa.o
	$(CC) $(CFLAGS) -o comparespa compare_nrel_spa.c $(OBJ) spa.o $(LFLAGS)

freespa.o: $(HDR) freespa.c	
clean:
	rm *.o testspa
