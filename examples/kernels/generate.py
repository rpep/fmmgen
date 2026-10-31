"""Generate C++ and Fortran 90 operators for a kernel that is not Laplace.

Run from this directory:  python generate.py
Writes yukawa_cpp.{h,cpp}, yukawa_f90.f90 and a JSON manifest for each.
"""
import fmmgen
from fmmgen.kernel import yukawa

kernel = yukawa(1.0)  # exp(-r)/r
orders = [2, 4, 6, 8]

for language, name in (("c++", "yukawa_cpp"), ("fortran", "yukawa_f90")):
    manifest = fmmgen.generate(kernel, max(orders), name, backend="cartesian",
                               orders=orders, language=language)
    print(name, manifest.files)
