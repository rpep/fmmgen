CXX=g++-9
CXXFLAGS=-Iinclude/ -fopenmp -Wall -fsanitize=address -fsanitize=undefined -fsanitize=leak -fno-omit-frame-pointer -fsignaling-nans -g
#CXXFLAGS=-Iinclude/ -fopenmp -Wall -O2 -std=c++14 -march=native

LIBNAME=fmm

SRC = src/operators.c   \
      src/calculate.cpp \
      src/tree.cpp      \
      src/utils.cpp


all: generate shared-library scaling reduced lazy direct


generate:
	python3 example.py

shared-library:
	${CXX} -c ${CXXFLAGS} -fpic src/operators.c   -o build/operators.o
	${CXX} -c ${CXXFLAGS} -fpic src/calculate.cpp -o build/calculate.o
	${CXX} -c ${CXXFLAGS} -fpic src/tree.cpp      -o build/tree.o
	${CXX} -c ${CXXFLAGS} -fpic src/utils.cpp     -o build/utils.o
	${CXX} -shared -o lib${LIBNAME}.so ${CXXFLAGS} build/operators.o \
                                                        build/calculate.o \
                                                        build/tree.o \
                                                        build/utils.o

static:
	${CXX} -c ${CXXFLAGS} src/operators.c -o build/operators.o
	${CXX} -c ${CXXFLAGS} src/calculate.cpp -o build/calculate.o
	${CXX} -c ${CXXFLAGS} src/tree.cpp -o build/tree.o
	${CXX} -c ${CXXFLAGS} src/utils.cpp -o build/utils.o

scaling-static: static
	${CXX} -c ${CXXFLAGS} src/apps/scaling_test.cpp scaling_test.o
	${CXX} ${CXXFLAGS} scaling_test.o build/*.o scaling

lazy-static: static
	${CXX} -c ${CXXFLAGS} src/apps/scaling_test_lazy.cpp -o lazy.o
	${CXX} ${CXXFLAGS} lazy.o build/*.o -o lazy-static


scaling: shared-library
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/scaling_test.cpp -o scaling -lfmm

reduced: shared-library
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/reduced_test.cpp -o reduced -lfmm

lazy: shared-library
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/scaling_test_lazy.cpp -o lazy -lfmm

direct: shared-library
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/direct.cpp -o direct -lfmm


matrix: shared-library
	${CXX} -Wl,-rpath,. -L. ${CXXFLAGS} src/apps/matrix.cpp -o matrix -lfmm

clean:
	rm -rf build/*
	rm -rf error*.txt
	rm -rf *.so
	rm -rf scaling
	rm -rf lib/*
	rm -rf *~
