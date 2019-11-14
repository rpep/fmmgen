# fmmgen
![Python version](https://img.shields.io/badge/Python-%3E%3D%203.6-brightgreen.svg)

This package generates Fast Multipole and Barnes-Hut operators for use in tree codes using automatic code generation.

It was written as part of the PhD research of Ryan Alexander Pepper at the
University of Southampton, under supervision of Hans Fangohr.

The library is written in Python (version 3.6). The package has few
dependencies; the main one is the SymPy library. Some parts of the SymPy library
are bundled within fmmgen due to changes needing to be made to the underlying
methods for the purposes of this code. Accordingly, fmmgen is licensed under the
3-Clause BSD License.

fmmgen consists of several parts:

1) Symbolical algebraic generation of the operators for fast multipole and Barnes-Hut codes in Cartesian Coordinates

Hand implementation of multipole formulae up to an arbitrary expansion order is non-trivial, and beyond 3rd order is a substantial effort. In general, this leaves most Cartesian fast multipole and Barnes-Hut authors writing operator functions by hand, as can be seen from the 

2) A code writer, which generates code from the expansion formulae. At present,
this generates C or C++ code but in future. The code makes use of something called 'common subexpression
elimination' (CSE) which reduces the number of operations which are performed in
the compiled code. Compilers already offer this functionality, but only at
higher levels of compiler optimisation is it turned on. Other optimisations are also present, for example, utilising the  

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

Once the package is installed, you can get started in using it to generate a C code version of the operators:

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

Alternatively, we suggest looking in the 'example' folder for a fully functioning OpenMP parallelised implementation of the FMM and Barnes-Hut methods using the code generated operators, which works for Coulomb, Dipole and higher order sources; all that needs to be done is change the 'source_order' parameter. By making other changes in the example.py file, one can enable or disable optimisations, which affects the run time significantly for some compilers. In general, we do not recommend the use of the GNU compiler, as in testing we find that the performance of the methods are significantly worse than when compiled with the Intel compiler. This has a side effect; we find that the symbolic algebra optimisations have less of an effect on the performance with the Intel compiler, which can factor expressions more effectively to avoid repeated computations than the GNU compiler at high optimisation levels.

## How to cite

A paper describing the work is under preparation. For now, please cite as

RA Pepper and H Fangohr, fmmgen, http://github.com/rpep/fmmgen (2019)

## References

The code was developed with particular reference to the following academic papers.

[1] Visscher, P. B., & Apalkov, D. M. (2010). Simple recursive implementation of fast multipole method. Journal of Magnetism and Magnetic Materials, 322(2), 275–281. https://doi.org/10.1016/j.jmmm.2009.09.033

[2] Coles, J. P., & Masella, M. (2015). The fast multipole method and point dipole moment polarizable force fields. The Journal of Chemical Physics, 142(2), 24109. https://doi.org/10.1063/1.4904922

[3] Beatson, R. and Greengard, L. (1997) A short course on fast multipole methods. In "Wavelets, Multilevel Methods and Elliptic PDEs", Oxford University Press, ISBN 0 19 850190 0

In addition, I would also point anyone interested in the Fast Multipole Method to the [video tutorial series](https://www.youtube.com/playlist?list=PLpa6_YduENMF080NikNninGG-7e1hK1eQ) of Dr. Rio Yokota of the Tokyo Institute of Technology for an overview and short course developing the 2-D method in a step-by-step way.

## Acknowledgements

This work was supported by the EPSRC Centre for Doctoral Training in Next
Generation Computational Modelling (grant EP/L015382/1).

I would also like to thank J. P. Coles for useful discussions regarding the method and implementation.
