#compiler and  option
FC=mpiifort
OPTF=-C -g -nofor-main

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CC = mpiicc
CFLAGS = -O3 -qopenmp -DAdd_

mumps=$(HOME)/Install/MUMPS_5.3.4
pord = $(mumps)/PORD
metis=$(HOME)/Install/metis-5.1.0/metis

LMKL= -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64

IPORD = -I$(pord)/include/
LPORD = -L$(pord)/lib/ -lpord
LMETIS = -L$(metis)/lib -lmetis
IMETIS = -I$(metis)/include
LMUMPS = -L$(mumps)/lib -lcmumps -ldmumps -lmumps_common
IMUMPS = -I$(mumps)/include


#===========================================================
BIN = .

LIB=$(LMUMPS) $(LMETIS) $(LMKL) $(LPORD)  -lm -qopenmp -pthread
INC=$(IMUMPS) $(IMETIS) $(IPORD)
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: clean main

##it has to be linked using mpiifort because main routines of MUMPS are in Fortran.
main:	$(OBJ)
	$(FC) $(OPTF) -o $(BIN)/main $(OBJ) $(LIB)

%.o: %.c 
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) main
