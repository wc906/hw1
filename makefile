cc=mpicc
flags=-lrt -lm


all: jacobi-mpi int_ring

jacobi-mpi: jacobi-mpi.c
	$(cc) $(flags) jacobi-mpi.c -o jacobi-mpi

int_ring: int_ring.c
	$(cc) $(flags) int_ring.c -o int_ring

array_ring: array_ring.c
	$(cc) $(flags) array_ring.c -o array_ring 

clean:
	rm jacobi-mpi int_ring array_ring