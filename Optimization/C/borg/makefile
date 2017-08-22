CC = gcc
ARGS = -Wall -O3
MPICC = mpicc
MPIARGS = -O3 -g

all:
	${CC} ${ARGS} -ansi -o mt19937ar.o -c mt19937ar.c
	${CC} ${ARGS} -ansi -o borg.o -c borg.c
	ar rcs libborg.a borg.o mt19937ar.o
	${CC} ${ARGS} -L. -o borg.exe frontend.c -lborg -lm
	${CC} ${ARGS} -ansi -L. -o dtlz2.exe dtlz2.c moeaframework.c -lm
	${CC} ${ARGS} -ansi -L. -o dtlz2_embedded.exe dtlz2_embedded.c -lborg -lm

mpi: all
	${MPICC} ${MPIARGS} -o borgms.o -c borgms.c
	ar rcs libborgms.a borgms.o mt19937ar.o
	${MPICC} ${MPIARGS} -L. -o dtlz2_ms.exe dtlz2_ms.c -lborgms -lm

clean: 
	rm -f mt19937ar.o borg.o libborg.a moeaframework.o borg.exe dtlz2.exe dtlz2_embedded.exe borgms.o libborgms.a dtlz2_ms.exe
