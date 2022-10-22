CC = gcc
CFLAGS = -O3 -g -Wall 

sxamg_lib = ../sxamg-master/src
sxamg_inc = ../sxamg-master/include


BIN = .
LIB = -lm -L$(sxamg_lib) -lsxamg
INC = -I. -I$(sxamg_inc)
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)



all: clean main

main:	$(OBJ)
	$(CC) $(CFLAGS) -o $(BIN)/main $(OBJ) $(LIB)

%.o: %.c 
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

clean:
	find . -name "*.o"   -exec rm {} \;
	find . -name "*.c%"  -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
	find . -name "*~"    -exec rm {} \;
	find . -name "\#*"   -exec rm {} \;
	rm -f $(OBJ) main


