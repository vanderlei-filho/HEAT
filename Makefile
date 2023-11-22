CC=mpicc
LD=mpicc

MPIDIR=${HOME}/opt/mpi
MPIINC=-I$(MPIDIR)/include
MPILIB=-lpthread -L$(MPIDIR)/lib -lmpi

# change this to point to your SCR install directory
SCRDIR=${HOME}/Downloads/scr-v3.0.1

SCRLIBDIR=-L$(SCRDIR)/install/lib -Wl,-rpath,$(SCRDIR)/install/lib -lscr
SCRINCLUDES=-I$(SCRDIR)/install/include

CFLAGS=-g -Wall
LDFLAGS= $(MPILIB) -g

LINK=$(LD)

APPS=jacobi_noft jacobi_ulfm jacobi_scr

all: $(APPS)

jacobi_noft: jacobi_noft.o main.o
	$(LINK) -o $@ $^ -lm

jacobi_ulfm: jacobi_ulfm.o main.o
	$(LINK) -o $@ $^ -lm

jacobi_scr: main_scr.o
	$(LINK) $(SCRINCLUDES) -o $@ $@.c $^ -lm $(SCRLIBDIR)

%.o: %.c jacobi.h
	$(CC) -c $(CFLAGS) -o $@ $<

clean:
	rm -f *.o $(APPS) *~
