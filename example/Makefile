
# Intel
CXX=icpc
CXXFLAGS=-O3 -std=c++14 -fopenmp -march=native
# Dangerous optimisations follow - here be dragons!
# CXXFLAGS=${CXXFLAGS}

# GCC
# CXX=g++
# Flags for debugging with GCC
# CXXFLAGS=-fsanitize=address -fno-omit-frame-pointer -g


all:
	rm -rf *.txt
	${CXX} ${CXXFLAGS} *.cpp -o main


