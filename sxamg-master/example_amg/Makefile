CC = gcc
CFLAGS = -O3 -g -Wall

sxamg_lib = ../src
sxamg_inc = ../include

BIN = .
LIB = -L$(sxamg_lib) -lsxamg -lm
INC = -I$(sxamg_inc) -I.
SRC = $(wildcard *.c)
OBJ = $(SRC:.c=.o)

all: clean main

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(INC) $(LIB)

main:	$(OBJ)
	$(CC) $(CFLAGS) -o $(BIN)/main $(OBJ) $(LIB)

clean:
	rm -f amg *.o
