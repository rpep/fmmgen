# Free-form Fortran 90. Any F90 compiler works; OpenMP is optional.
FC=gfortran
FFLAGS=-O3 -fopenmp
# Strict-standard check:  make FFLAGS="-O2 -std=f2003 -Wall -Wextra"

all: main

operators.f90: example.py
	uv run python example.py

operators.o: operators.f90
	${FC} ${FFLAGS} -c operators.f90

main: operators.o fmm_example.f90
	${FC} ${FFLAGS} operators.o fmm_example.f90 -o main

clean:
	rm -f *.o *.mod main
