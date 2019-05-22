CXX=g++-8
CXXFLAGS=-Iinclude/ -fopenmp -Wall

LIBNAME=bh

SRC = src/operators.c   \
      src/calculate.cpp \
      src/tree.cpp      \
      src/utils.cpp


all: generate library scaling


generate:
	# python3 example.py

library:
	${CXX} -c ${CXXFLAGS} -fpic src/operators.c   -o build/operators.o
	${CXX} -c ${CXXFLAGS} -fpic src/calculate.cpp -o build/calculate.o
	${CXX} -c ${CXXFLAGS} -fpic src/tree.cpp      -o build/tree.o
	${CXX} -c ${CXXFLAGS} -fpic src/utils.cpp     -o build/utils.o
	${CXX} -shared -o lib${LIBNAME}.so ${CXXFLAGS} build/operators.o \
                                                        build/calculate.o \
                                                        build/tree.o \
                                                        build/utils.o

scaling: library
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/scaling_test.cpp -o scaling -lbh
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/dipole_field.cpp -o dipole -lbh


clean:
	rm -rf build/*
	rm -rf error*.txt
	rm -rf *.so
	rm -rf scaling
	rm -rf dipole
	rm -rf lib/*
	rm -rf *~