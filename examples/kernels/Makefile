# Generic-kernel example: generate operators for exp(-r)/r, then build and run
# the C++ and Fortran 90 checks. Needs g++; the Fortran target also needs gfortran.
PYTHON?=python
CXX=g++
CXXFLAGS=-O2 -std=c++14
FC=gfortran
FFLAGS=-O2

all: run_cpp

yukawa_cpp.cpp yukawa_cpp.h yukawa_f90.f90: generate.py
	${PYTHON} generate.py

single_cell_cpp: single_cell.cpp yukawa_cpp.cpp yukawa_cpp.h
	${CXX} ${CXXFLAGS} single_cell.cpp yukawa_cpp.cpp -o $@

single_cell_f90: single_cell.f90 yukawa_f90.f90
	${FC} ${FFLAGS} yukawa_f90.f90 single_cell.f90 -o $@

run_cpp: single_cell_cpp
	./single_cell_cpp

run_f90: single_cell_f90
	./single_cell_f90

clean:
	rm -f yukawa_cpp.* yukawa_f90.* *.mod single_cell_cpp single_cell_f90
