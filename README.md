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

The code writer can also output a Cython wrapper for this C code, which can be
used for quick testing of the operators.


To try out the module:

```python
import fmmgen


order = 4
module_name = "operators"

# To get the unoptimized version of the code, generated without CSE:
fmmgen.generate_code(order, module_name, generate_cython_wrapper=True, CSE=False) 

# Alternatively, for the optimised version:
fmmgen.generate_code(order, module_name, generate_cython_wrapper=True, CSE=True)

# When Cython generation is enabled, it is possible to use the operator functions:

import pyximport
pyximport.install()
import operators_wrap

```