
# Intel
# CXX=icpc
# CXXFLAGS=-O3 -std=c++14 -fopenmp -march=native
# Dangerous optimisations follow - here be dragons!
# CXXFLAGS=${CXXFLAGS}

# GCC from Homebrew
CXX=g++-15
CXXFLAGS=-O3 -std=c++14 -fopenmp -march=native
# Flags for debugging with GCC
# CXXFLAGS=-fsanitize=address -fno-omit-frame-pointer -g


# The generated operators are long straight-line scalar code. Measured here, -O1 compiles
# them in about 2/3 of the -O3 time and runs the same (M2L is scalar FMA throughput either
# way), so they get their own, cheaper flags. Set OPFLAGS=${CXXFLAGS} to compile them like
# everything else.
OPFLAGS=-O1 -std=c++14 -fopenmp -march=native

PYTHON?=python
SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)

# `make -j` compiles the translation units (and every operators_<operator>.cpp) in parallel.
all: main

main: $(OBJ)
	rm -rf *.txt
	${CXX} -fopenmp $(OBJ) -o main

operators.h:
	@echo "operators.h is missing: generate the operators first with 'make generate' (or python example.py)"; false

OPOBJ=$(filter operators%.o,$(OBJ))
$(OPOBJ): %.o: %.cpp operators.h
	${CXX} ${OPFLAGS} -c $< -o $@

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c $< -o $@

# Regenerate the operators using every CPU for SymPy and write each operator to its own file so that
# `make -j` can compile them in parallel. FMMGEN_CACHE reuses already-generated orders.
generate:
	FMMGEN_JOBS=0 FMMGEN_SPLIT=1 FMMGEN_CACHE=.fmmgen_cache ${PYTHON} example.py

clean:
	rm -f *.o main


