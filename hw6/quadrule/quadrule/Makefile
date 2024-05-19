#-----------------------------------------------------------

# NONLINEAR_SOLVER={NEWTON,MINPACK,PETSC,TENSOLVE,CCMNF,GSLSIMAN}
CDEFS0=-DNONLINEAR_SOLVER=MINPACK	
CDEFS = -DUSE_MPI=0 ${CDEFS0}

#MPFR_CPPFLAGS = -DUSE_MPFR=1
#MPFR_LIB = -lmpfr -lgmp

MINPACK_LIB = minpack/libminpack.a
GSL_CPPFLAGS =
GSL_LIB = ##-lgsl

GAUL_CPPFLAGS =
GAUL_LIB = #-lgaul -lgaul_util

# PETSc
# Headers (-I${PETSC_DIR}/include -I${PETSC_DIR}/${PETSC_ARCH}/include)
PETSC_CPPFLAGS = $(shell sh petsc-configs.sh CPPFLAGS)
# Libs (-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -llapack -lblas -lpthread)
PETSC_LIBS = $(shell sh petsc-configs.sh LDFLAGS RPATHS LIBS)

CC = gcc
CFLAGS = -Wall -O2 -g
CPPFLAGS = ${PETSC_CPPFLAGS}
F77 = gfortran
# Note: with GCC 10 '-std=legacy' is needed, otherwise 'Rank mismatch in argument' errors are produced when compiling the original toms-768/Src/Dp/src.f file.
FFLAGS = -g -O2 #-std=legacy
LDFLAGS =
LINKER = gfortran
LIBS = ${PETSC_LIBS} #-lquadmath

quadrule: quadrule.o
	@(cd minpack; ${MAKE} F77="$(F77)" FFLAGS="$(FFLAGS)" libminpack.a)
	@(cd ccmnf; ${MAKE} CC="$(CC)" CFLAGS="$(CFLAGS)" LINKER="$(LINKER)" ccmNF.o)
	${LINKER} -o $@ $^ ccmnf/ccmNF.o ${LDFLAGS} ${MINPACK_LIB} ${GAUL_LIB} ${GSL_LIB} ${MPFR_LIB} ${LIBS}

quadrule.o: quadrule.c quadrule.h Makefile
	@if ! test -r ccmnf; then ln -s ../misc/ccmnf .; fi
	${CC} ${CFLAGS} ${CPPFLAGS} ${GAUL_CPPFLAGS} ${GSL_CPPFLAGS} ${MPFR_CPPFLAGS} ${CDEFS} -c quadrule.c

clean:
	rm -f *.o quadrule {minpack,ccmnf}/*.[oa]
