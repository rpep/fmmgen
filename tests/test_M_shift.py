import fmmgen.generator as gen
import fmmgen.expansions as exp
from fmmgen.utils import Nterms
import sympy as sp

x, y, z = sp.symbols("x y z")
symbols = (x, y, z)


def test_M_shift_monopole():
    order = 2
    source = 0

    M = sp.MatrixSymbol("M", Nterms(order), 1)
    M_dict, _ = gen.generate_mappings(
        order, symbols, key="grevlex", source_order=source
    )
    # Check monopole term
    assert exp.M_shift((0, 0, 0), order, symbols, M_dict, source_order=source) == M[0]
    # Check dipole terms
    assert (
        exp.M_shift((1, 0, 0), order, symbols, M_dict, source_order=source)
        == M[1] + M[0] * x
    )
    assert (
        exp.M_shift((0, 1, 0), order, symbols, M_dict, source_order=source)
        == M[2] + M[0] * y
    )
    assert (
        exp.M_shift((0, 0, 1), order, symbols, M_dict, source_order=source)
        == M[3] + M[0] * z
    )
    # Check quadrupole terms
    assert (
        exp.M_shift((2, 0, 0), order, symbols, M_dict, source_order=source)
        == x**2 / 2 * M[0] + x * M[1] + M[4]
    )

    assert (
        sp.expand(exp.M_shift((1, 1, 0), order, symbols, M_dict, source_order=source))
        == x * y * M[0] + y * M[1] + x * M[2] + M[5]
    )

    assert (
        sp.expand(exp.M_shift((1, 0, 1), order, symbols, M_dict, source_order=source))
        == x * z * M[0] + z * M[1] + x * M[3] + M[6]
    )


def test_M_shift_dipole():
    order = 2
    source = 1
    M = sp.MatrixSymbol("M", Nterms(order), 1)
    M_dict, _ = gen.generate_mappings(
        order, symbols, key="grevlex", source_order=source
    )

    # With a source order of 1, the monopole term is no longer produced. We should check
    # therefore, that M_shift returns an error.

    print(M_dict[(1, 0, 1)])

    try:
        exp.M_shift((0, 0, 0), order, symbols, M_dict, source_order=source)
        raise AssertionError("No exception raised!")
    except ValueError:
        pass

    assert exp.M_shift((1, 0, 0), order, symbols, M_dict, source_order=source) == M[0]
    assert exp.M_shift((0, 1, 0), order, symbols, M_dict, source_order=source) == M[1]
    assert exp.M_shift((0, 0, 1), order, symbols, M_dict, source_order=source) == M[2]

    assert (
        sp.expand(exp.M_shift((2, 0, 0), order, symbols, M_dict, source_order=source))
        == M[3] + M[0] * x
    )

    assert (
        sp.expand(exp.M_shift((1, 1, 0), order, symbols, M_dict, source_order=source))
        == y * M[0] + x * M[1] + M[4]
    )

    assert (
        sp.expand(exp.M_shift((1, 0, 1), order, symbols, M_dict, source_order=source))
        == z * M[0] + x * M[2] + M[5]
    )


def test_generate_M_shift_operators_monopole():
    order = 1

    M_dict = {(0, 0, 0): 0, (1, 0, 0): 1, (0, 1, 0): 2, (0, 0, 1): 3}

    M = sp.MatrixSymbol("M", 4, 1)
    x, y, z = sp.symbols("x y z")

    answer = [M[0], M[1] + M[0] * x, M[2] + M[0] * y, M[3] + M[0] * z]
    result = gen.generate_M_shift_operators(order, (x, y, z), M_dict)
    for a, b in zip(answer, result):
        assert a == b


def test_generate_M_shift_operators_dipole():
    order = 2
    source = 1
    M_dict = {
        (1, 0, 0): 0,
        (0, 1, 0): 1,
        (0, 0, 1): 2,
        (2, 0, 0): 3,
        (1, 1, 0): 4,
        (1, 0, 1): 5,
        (0, 2, 0): 6,
        (0, 1, 1): 7,
        (0, 0, 2): 8,
    }

    M = sp.MatrixSymbol("M", Nterms(order), 1)
    x, y, z = sp.symbols("x y z")

    # Hand computed for correctness
    answer = [
        M[0],
        M[1],
        M[2],
        M[0] * x + M[3],
        M[0] * y + M[1] * x + M[4],
        M[0] * z + M[2] * x + M[5],
        M[1] * y + M[6],
        M[1] * z + M[2] * y + M[7],
        M[2] * z + M[8],
    ]

    result = gen.generate_M_shift_operators(
        order, (x, y, z), M_dict, source_order=source
    )
    for i, (a, b) in enumerate(zip(answer, result)):
        assert a == b, f"i = {i}, a = {a}, b = {b}"
