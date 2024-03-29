# fmmgen
![Python version](https://img.shields.io/badge/Python-%3E%3D%203.10-brightgreen.svg)
![C++14 version](https://img.shields.io/badge/c%2B%2B-14-brightgreen)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3842591.svg)](https://doi.org/10.5281/zenodo.3842591)
[![Arxiv Paper](https://img.shields.io/badge/arxiv-2005.12351-B31B1B)](https://arxiv.org/abs/2005.12351)
[![.github/workflows/actions.yml](https://github.com/rpep/fmmgen/actions/workflows/actions.yml/badge.svg)](https://github.com/rpep/fmmgen/actions/workflows/actions.yml)

This package generates Fast Multipole and Barnes-Hut operators for use in tree codes.
It was written as part of the PhD research of Ryan Alexander Pepper at the University of Southampton.

The library is written in Python, and requires at least version 3.6. The package has few dependencies; the main one is the SymPy library. Some parts of the SymPy library are bundled within fmmgen due to changes needing to be made to the underlying methods for the purposes of this code. Accordingly, fmmgen is licensed under the 3-Clause BSD License.

fmmgen consists of several parts:

1) Symbolical algebraic generation of the operators for fast multipole and Barnes-Hut codes in Cartesian Coordinates. Hand implementation of multipole formulae up to an arbitrary expansion order is non-trivial, and beyond 3rd order is a substantial effort. In general, this leaves most Cartesian fast multipole and Barnes-Hut authors writing operator functions by hand, as can be seen from other similar packages.

2) A code writer, which generates code from the expansion formulae. At present,
this generates C or C++ code but in future may be extended. The code makes use of something called 'common subexpression
elimination' (CSE) which reduces the number of operations which are performed in
the compiled code. Compilers already offer this functionality, but only at
higher levels of compiler optimisation is it turned on. Other optimisations are also present.

The code writer can also output a Cython wrapper for this C or C++ code, which can be
used for quick testing of the operators.


## Installation

To try out the module, first install it and the requirements:

```bash
git clone https://github.com/rpep/fmmgen.git
cd fmmgen
pip install .
```

## Example

Once the package is installed, you can get started in using it to generate a C code version of the operators:

```python
import fmmgen
import numpy as np

# Order of the multipole expansion
order = 10
# Order of the sources (i.e. monopole = 0, dipole = 1)
source_order = 0

# This module name is used to label the source files, so
# with this, we get output of 'operators.c' and 'operators.h'
module_name = "operators"

# Set whether ultimately, the field or potential are to be calculated.
# Note that calculating the field constrains the operator functions;
# no 0th order function is generated for the particle-to-multipole
# (P2M) operator.
potential = True
field = False

# Choose a language ('c' or 'cpp')
language = 'c'

# Choose whether to enable the common-subexpression elimination optimisation:
CSE = True

# Choose at what order expressions like x*x*x*x*x are converted to pow(x, 5)
minpow = 5

# Enable/disable generation of Cython code to allow the use of the functions from Python:
cython=True

# To get the unoptimized version of the code, generated without common subexpression elimination:
fmmgen.generate_code(order, module_name, cython=cython, CSE=CSE,
                     source_order=source_order,
                     minpow=minpow,
                     potential=True,
                     field=False)


# When Cython generation is enabled, it is possible to use the operator functions
# directly from Python by importing them with pyximport. These have the same
# API as the C code that is exported.
import pyximport
pyximport.install()
import operators_wrap as fmm

# To calculate the multipole moments of a charge q located at (0, 0, d)
# about the origin, for example, you can use the following:
d = 2.0
q = 3.0
# Position of the charge:
r = np.array([0.0, 0.0, d])
# Number of entries in a multipole array for quadrupoles:
Nterms = fmmgen.utils.Nterms(2)

# Multipole input array:
Q = np.zeros(Nterms)
Q[0] = q
# Multipole output array:
M = np.zeros(Nterms)
fmm.M2M(*r, Q, M, 2)
print(M)
# Expected output:
#
# M = [3.0,            [monopole moment]
#      0.0, 0.0, 6.0,  [x, y, z dipole moments]
#      0.0, 0.0, 0.0,
#      0.0, 0.0, 6.0]  [xx, xy, xz, yy, yz, zz quadrupole moments]
#
```

We suggest looking in the 'example' folder for a fully functioning OpenMP parallelised implementation of the FMM and Barnes-Hut methods using the code generated operators, which works for Coulomb, Dipole and higher order sources; all that needs to be done is change the 'source_order' parameter. By making other changes in the example.py file, one can enable or disable optimisations, which affects the run time significantly for some compilers.

In general, we do not recommend the use of the GNU compiler for this; in testing we find that the performance of the methods are significantly worse than when compiled with the Intel compiler. This has a side effect; we find that the symbolic algebra optimisations have less of an effect on the performance with the Intel compiler, which can factor expressions more effectively to avoid repeated computations than the GNU compiler at high optimisation levels.

## References

The code was developed with particular reference to the following academic papers.

[1] Visscher, P. B., & Apalkov, D. M. (2010). Simple recursive implementation of fast multipole method. Journal of Magnetism and Magnetic Materials, 322(2), 275–281. https://doi.org/10.1016/j.jmmm.2009.09.033

[2] Coles, J. P., & Masella, M. (2015). The fast multipole method and point dipole moment polarizable force fields. The Journal of Chemical Physics, 142(2), 24109. https://doi.org/10.1063/1.4904922

[3] Beatson, R. and Greengard, L. (1997) A short course on fast multipole methods. In "Wavelets, Multilevel Methods and Elliptic PDEs", Oxford University Press, ISBN 0 19 850190 0

In addition, I would also point anyone interested in the Fast Multipole Method to the [video tutorial series](https://www.youtube.com/playlist?list=PLpa6_YduENMF080NikNninGG-7e1hK1eQ) of Dr. Rio Yokota of the Tokyo Institute of Technology for an overview and short course developing the 2-D method in a step-by-step way.

I would also like to thank J. P. Coles for useful discussions regarding the method and implementation.

## Citations

The following papers have cited or used Fmmgen:

[1] [Efficient Open-Source Implementations of Linear-Scaling Polarizable Embedding: Use Octrees to Save the Trees](https://doi.org/10.1021/acs.jctc.1c00225) M. Scheurer, P. Reinholdt, J. M. H. Olsen, A. Dreuw, J Kongsted, J. Chem. Theory Comput. 17, 6, 3445–3454 (2021)
