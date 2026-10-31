"""The public entry points: bad input is rejected with a clear error, and generate_code(kernel=...) works."""
import pytest
import sympy as sp

import fmmgen
from fmmgen.kernel import Kernel, generate_kernel_code, power_law, yukawa

x, y, z = sp.symbols("x y z")
R = sp.sqrt(x**2 + y**2 + z**2)


def test_kernel_rejects_bad_definitions():
    a, b, c, d = sp.symbols("a b c d")
    with pytest.raises(ValueError, match="dimension"):
        Kernel(1 / sp.sqrt(a**2 + b**2 + c**2 + d**2), (a, b, c, d))
    with pytest.raises(ValueError, match="unbound symbols"):
        Kernel(sp.exp(-sp.Symbol("k") * R) / R, (x, y, z))
    with pytest.raises(ValueError, match="method must be"):
        Kernel(1 / R, (x, y, z), method="bogus")
    with pytest.raises(ValueError, match="need a radial kernel"):
        Kernel(1 / (R + x), (x, y, z), method="recurrence")
    # a pde the kernel does not satisfy would silently give wrong derivatives: it is checked
    with pytest.raises(ValueError, match="does not satisfy"):
        Kernel(1 / R**3, (x, y, z), pde={(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 2): 1})
    with pytest.raises(ValueError, match="multi-index length"):
        Kernel(1 / R, (x, y, z), pde={(2, 0): 1})
    yukawa(0.7)  # and a correct one is accepted


def test_generate_kernel_code_rejects_bad_options(tmp_path):
    kw = dict(src_dir=str(tmp_path), include_dir=str(tmp_path))
    with pytest.raises(ValueError, match="language"):
        generate_kernel_code(power_law(1), 2, "k", language="cobol", **kw)
    with pytest.raises(ValueError, match="source_order"):
        generate_kernel_code(power_law(1), 2, "k", orders=[1, 2], source_order=2, **kw)


def test_generate_code_delegates_to_the_kernel_generator(tmp_path):
    fmmgen.generate_code(3, "yk", kernel=yukawa(0.7), language="c", src_dir=str(tmp_path),
                         include_dir=str(tmp_path))
    header = (tmp_path / "yk.h").read_text()
    assert "M2L_3" in header and "FMMGEN_MAXORDER 3" in header
    # minpow=0 (generate_code's default) means automatic for a kernel: no pow() calls in the operators
    assert "pow(" not in (tmp_path / "yk.c").read_text()
