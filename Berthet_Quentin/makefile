CC=mpicc
CC_FLAGS=-std=gnu11 -Wall -Wextra -Wpedantic -std=c99 -MMD -g -fsanitize=undefined
L=-lm
NUM_PROC=4
MPI=mpirun -np $(NUM_PROC) ./mmm.out

all: clean build

run: clean build
	$(MPI)

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) $(L)

build: mmm.o para_mmm.o mmm_utils.o 
	$(CC) -o mmm.out mmm.o para_mmm.o mmm_utils.o $(CFLAGS) $(L)


clean:
	rm -f mmm.out mmm.o para_mmm.o mmm_utils.o



