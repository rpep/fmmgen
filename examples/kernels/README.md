# Generic-kernel example

Operators for a kernel other than Laplace, generated through `fmmgen.generate` and checked against direct summation from C++ and from Fortran 90.

The kernel is the Yukawa potential exp(-r)/r. Change `kernel` in `generate.py` to any `fmmgen.kernel` library kernel or your own `Kernel`.

```
PYTHONPATH=../.. make run_cpp     # needs g++
PYTHONPATH=../.. make run_f90     # needs gfortran
```

`PYTHON=<path>` selects the interpreter that has sympy and numpy installed.

## What it does

`single_cell.cpp` and `single_cell.f90` build two source cells and two target cells around a parent cell each, then run the full chain S2M, M2M, M2L, L2L, L2P and compare the potential with the generated P2P. They use the same random number generator, so the two programs see the same points and print the same errors.

Expected output (the error falls with the order):

```
 order  max relative error in potential
     2  2.159e-01
     4  9.608e-03
     6  1.871e-04
     8  1.079e-05
```

## Conventions

All operators accumulate with `+=`. Displacements are S2M(centre−source), M2M(parent−child), M2L(target−source), L2L(child−parent), L2P(body−centre) and P2P(target−source). The output is {potential, −∇potential}. The size of each array is in `yukawa_*.json` (`layouts`).

The `{name}.json` manifest lists the files, functions and array sizes, so a driver can check what was generated.
