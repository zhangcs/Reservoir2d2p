
.SUFFIXES: .c .f90 .f .for

FC = gfortran
CC = gcc

C_COMPILE_FLAGS = -g -Wall # -O2
F_COMPILE_FLAGS = -g -Wall # -O2

CFLAGS = ${C_COMPILE_FLAGS}
FFLAGS = ${F_COMPILE_FLAGS}

.c.o:
	${CC} -o $@ -c ${CFLAGS} $<
.f.o:
	${FC} -o $@ -c ${FFLAGS} $<
.for.o:
	${FC} -o $@ -c ${FFLAGS} $<
.f90.o:
	${FC} -o $@ -c ${FFLAGS} $<

OBJS =  main.o \
	oilgas.o \
	iniguess.o \
	newton.o \
	pdterms.o \
	sdterms.o \
	trans.o \
	residual.o \
	jacobi.o \
	solver.o \
	tool.o \

test :  ${OBJS}
	${CC} -o 2d2p ${OBJS} -lm 

clean :
	-rm -f *.o *~

cleanout:
	-rm -rf 2d2p ./output/*

allclean:
	make clean
	make cleanout

