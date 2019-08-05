# fmmgen
![Python version](https://img.shields.io/badge/Python-%3E%3D%203.6-brightgreen.svg)

This package generates Fast Multipole and Barnes-Hut operators for use in tree codes.
It was written as part of the PhD research of Ryan Pepper.

It consists of several parts:

1) Symbolical algebraic generation of the operators.
Hand implementation of multipole formulae up to an arbitrary expansion order is
non-trivial, and beyond 3rd order is a substantial effort.

2) A code writer, which generates code from the expansion formulae. At present,
this generates C code but in future, the option to generate Fortran functions
will be added. The code makes use of something called 'common subexpression
elimination' (CSE) which reduces the number of operations which are performed in
the compiled code. Compilers already offer this functionality, but only at
higher levels of compiler optimisation is it turned on.

The code writer can also output a Cython wrapper for this C or C++ code, which can be
used for quick testing of the operators.


## Installation

To try out the module, first install it and the requirements:

```bash
git clone https://github.com/rpep/fmmgen.git
cd fmmgen
pip install -r requirements.txt
pip install .
```

## Example

Then you can get started in using it to generate a C code version of the operators:

```python
import fmmgen

order = 4
module_name = "operators"

# To get the unoptimized version of the code, generated without CSE:
fmmgen.generate_code(order, module_name, cython=True, CSE=False)

# Alternatively, for the optimised version:
fmmgen.generate_code(order, module_name, cython=True, CSE=True)

# When Cython generation is enabled, it is possible to use the operator functions
# by importing them with pyximport:
import pyximport
pyximport.install()
import operators_wrap as fmm
```

Alternatively, after generating the code as above, the generated C code can be compiled and used directly (following instructions should work on Linux/MacOS):

```
// main.c
#include "stdio.h"
#include "stdlib.h"

#include "operators.h"

unsigned int Nterms(const unsigned int order) {
  int nterms = 0;
  for(unsigned int n = 0; n <= order; n++) {
    nterms += (n*(n + 1)) / 2;
  }
  return nterms;
}

int main() {
  // Expansion Order:
  const unsigned int order = 2;
  // Number of terms in expansion at order:
  unsigned int N = Nterms(order);
  // Single multipole array set to zero:

  double *M = (double *) calloc(N, sizeof(double));

  double x, y, z, q;
  x = 1.1;
  y = 2.4;
  z = -1.2;
  q = 1e-5;

  // Calculate P2M expansion up to order
  P2M(x, y, z, q, M, order);

  for(unsigned int i = 0; i < N; i++) {
    printf("M[%d] = %g\n", i, M[i]);
  }


  free(M);
  return 0;
}
```

To compile and run:
```
gcc -c operators.c
gcc -c main.c
gcc operators.o main.o -o main
./main
```

See the folder 'example' for a fully working, OpenMP Barnes-Hut and Fast Multipole Method code, built using the operators.

## Tests

To run the tests, simply run from the project root directory:

```
make runtests
```

## References

The code was developed with reference to the following academic papers.

[1] Visscher, P. B., & Apalkov, D. M. (2010). Simple recursive implementation of fast multipole method. Journal of Magnetism and Magnetic Materials, 322(2), 275–281. https://doi.org/10.1016/j.jmmm.2009.09.033

[2] Coles, J. P., & Masella, M. (2015). The fast multipole method and point dipole moment polarizable force fields. The Journal of Chemical Physics, 142(2), 24109. https://doi.org/10.1063/1.4904922

I would also like to thank J. P. Coles for useful discussions regarding the method and implementation.
