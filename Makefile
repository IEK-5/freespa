SRC=freespa.c testspa.c
OBJ=freespa.o
HDR=freespa.h freespa_tables.h freespa_dt_table.h
CC=gcc
CFLAGS=-O3 -flto
# CFLAGS=-Og -Wall -pedantic -flto -g
LFLAGS=-lm
ifneq ("$(wildcard spa.c)","")
	ifneq ("$(wildcard spa.h)","")
		NRELSPA = -DNRELSPA
		OBJ:=$(OBJ) spa.o
	endif
endif
$(info $$NRELSPA is [${NRELSPA}])

all: freespa.o check
distdir:
	cp $(SRC) $(HDR) Makefile reference.dat $(distdir)
check: test
	./testspa t reference.dat
	./testtime
test: $(SRC) $(HDR) $(OBJ) testspa.c testtime.c
	$(CC) $(CFLAGS) -o testspa testspa.c $(OBJ) $(LFLAGS)
	$(CC) $(CFLAGS) -o testtime testtime.c $(OBJ) $(LFLAGS)

compare: $(SRC) $(HDR) $(OBJ) spa.c spa.h
	$(CC) $(CFLAGS) -o comparespa compare_nrel_spa.c $(OBJ) $(LFLAGS)
	
spa: freespa.o cli_spa.c $(OBJ)
	$(CC) $(CFLAGS) -o spa $(NRELSPA) cli_spa.c $(OBJ) $(LFLAGS)

freespa.o: $(HDR) freespa.c	
clean:
	rm *.o testspa
# targets for autotools compatibility
install:
install-data:
install-exec:
uninstall:
install-dvi:
install-html:
install-info:
install-ps:
install-pdf:
installdirs:
installcheck:
mostlyclean:
distclean: clean
maintainer-clean: clean
dvi:
pdf:
ps:
info:
html:
tags:
ctags:
