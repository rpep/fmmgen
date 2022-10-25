import fmmgen.generator as gen
import fmmgen.expansions as exp
from fmmgen.utils import q
import sympy as sp

x, y, z = sp.symbols("x y z")
symbols = (x, y, z)


def test_M_monopole():
    order = 2
    source = 0

    M_dict, _ = gen.generate_mappings(
        order, symbols, key="grevlex", source_order=source
    )
    # Check monopole term
    assert exp.M((0, 0, 0), symbols) == q
    # Check dipole terms
    assert exp.M((1, 0, 0), symbols) == -q * x
    assert exp.M((0, 1, 0), symbols) == -q * y
    assert exp.M((0, 0, 1), symbols) == -q * z


def test_M_dipole():
    order = 2
    source = 1
    M_dict, _ = gen.generate_mappings(
        order, symbols, key="grevlex", source_order=source
    )

    mux, muy, muz = sp.symbols("mux muy muz")

    # Check dipole terms
    assert exp.M_dipole((1, 0, 0), symbols, M_dict) == mux
    assert exp.M_dipole((0, 1, 0), symbols, M_dict) == muy
    assert exp.M_dipole((0, 0, 1), symbols, M_dict) == muz

    # Check quadrupole terms
    assert exp.M_dipole((2, 0, 0), symbols, M_dict) == -mux * x
    assert exp.M_dipole((1, 1, 0), symbols, M_dict) == -mux * y - muy * x
    assert exp.M_dipole((1, 0, 1), symbols, M_dict) == -mux * z - muz * x
    assert exp.M_dipole((0, 2, 0), symbols, M_dict) == -muy * y
    assert exp.M_dipole((0, 1, 1), symbols, M_dict) == -muy * z - muz * y
    assert exp.M_dipole((0, 0, 2), symbols, M_dict) == -muz * z
