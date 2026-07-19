CC      = mpicc
CFLAGS  = -O2 -fopenmp -Wall -Wextra
LDFLAGS = -lm

NP      ?= 4
THREADS ?= 2

all: heat2d

heat2d: heat2d.c
	$(CC) $(CFLAGS) -o heat2d heat2d.c $(LDFLAGS)

run: heat2d
	OMP_NUM_THREADS=$(THREADS) mpirun -np $(NP) ./heat2d entrada.txt saida

seq: heat2d
	OMP_NUM_THREADS=1 mpirun -np 1 ./heat2d entrada.txt saida_seq

clean:
	rm -f heat2d
	rm -rf saida saida_seq

.PHONY: all run seq clean
