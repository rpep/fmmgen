CXX=g++-8
CXXFLAGS=-Iinclude/ -fopenmp -Wall -fsanitize=address -fsanitize=undefined -fsanitize=leak -fno-omit-frame-pointer -fsignaling-nans -g
#CXXFLAGS=-Iinclude/ -fopenmp -Wall -O3 -std=c++14

LIBNAME=fmm

SRC = src/operators.c   \
      src/calculate.cpp \
      src/tree.cpp      \
      src/utils.cpp


all: generate library scaling reduced lazy


generate:
#	python3 example.py

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
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/scaling_test.cpp -o scaling -lfmm

reduced: library
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/reduced_test.cpp -o reduced -lfmm

lazy: library
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/scaling_test_lazy.cpp -o lazy -lfmm

clean:
	rm -rf build/*
	rm -rf error*.txt
	rm -rf *.so
	rm -rf scaling
	rm -rf lib/*
	rm -rf *~
