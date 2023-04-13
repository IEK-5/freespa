SRC=freespa.c testspa.c
OBJ=freespa.o
HDR=freespa.h freespa_tables.h freespa_dt_table.h
CC=gcc
CFLAGS=-O3 -flto
# my tests with Ofast suggest no problems, shaves about 40% off the execution time
#CFLAGS=-Ofast -flto -march=native
# CFLAGS=-Og -Wall -pedantic -flto -g
# CC=x86_64-w64-mingw32-gcc
# CFLAGS=-O3 -flto -D_POSIX_C_SOURCE -D_POSIX_THREAD_SAFE_FUNCTIONS
# CFLAGS=-Og -Wall -pedantic -flto -g -D_POSIX_C_SOURCE
LFLAGS=-lm
ifneq ("$(wildcard spa.c)","")
	ifneq ("$(wildcard spa.h)","")
		NRELSPA = -DNRELSPA
		OBJ:=$(OBJ) spa.o
	endif
endif
deflat=50.90329
deflon=6.41143
defele=96.0


all: freespa.o check spa
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
	$(CC) $(CFLAGS) -o spa $(NRELSPA) -DDEFLAT=$(deflat) -DDEFLON=$(deflon) -DDEFELE=$(defele) cli_spa.c $(OBJ) $(LFLAGS)

freespa.o: $(HDR) freespa.c	
clean:
	$(RM) *.o testspa testtime spa comparespa
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
