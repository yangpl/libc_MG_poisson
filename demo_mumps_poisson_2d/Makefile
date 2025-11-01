# GNU Compiler and options, assume metis,scotch,lapack,scalapack have been install using sudo apt-get install libmetis-dev libparmetis-dev libscotch-dev libptscotch-dev libatlas-base-dev openmpi-bin libopenmpi-dev liblapack-dev libscalapack-openmpi-dev
FC = mpif90 
CC = mpicc 

CFLAGS = -O3 -fopenmp -DAdd_ #validate Add_=1 in the code
FFLAGS = -g -fopenmp

#===========================================================
# MUMPS paths
mumps = $(HOME)/Install/MUMPS_5.7.3
IMUMPS = -I$(mumps)/include
LMUMPS = -L$(mumps)/lib -ldmumps -lmumps_common

pord = $(mumps)/PORD
IPORD = -I$(pord)/include
LPORD = -L$(pord)/lib -lpord

IMETIS    = -I/usr/include/parmetis
LMETIS = -L/usr/lib  -lparmetis -lmetis #parallel version
#IMETIS    = -I/usr/include/metis #sequential version

ISCOTCH   = -I/usr/include/scotch
LSCOTCH   = -L/usr/lib -lptesmumps -lptscotch -lptscotcherr #parallel version
#LSCOTCH   = -L/usr/lib -lesmumps -lscotch -lscotcherr #sequential version

#===========================================================
BIN = .

# Order: MUMPS -> PORD -> SCOTCH -> ScaLAPACK -> LAPACK -> BLAS -> system libraries
LIB = $(LMUMPS) $(LPORD) $(LMETIS) $(LSCOTCH) -lscalapack-openmpi -llapack -lblas -lpthread -lm 
INC = $(IMUMPS) $(IPORD) $(LMETIS) $(ISCOTCH)

SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: clean main

main: $(OBJ)
	$(FC) $(FFLAGS) -o $(BIN)/main $(OBJ) $(LIB)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) main
