import fmmgen
import numpy as np
from fmmgen.utils import Nterms
import logging

logging.getLogger("fmmgen")

TOTALORDER = 4


def test_linear_dipole():
    source_order = 0
    order = source_order + TOTALORDER
    cse = True
    atomic = True
    precision = "double"

    fmmgen.generate_code(
        order,
        "LinearDipole",
        precision=precision,
        CSE=cse,
        cython=True,
        potential=True,
        field=True,
        source_order=source_order,
        atomic=atomic,
        minpow=5,
        harmonic_derivs=True,
        language="c",
    )

    import pyximport

    pyximport.install()
    import LinearDipole_wrap as fmm

    # Locate charge +q at (0, 0, d/2)
    # Locate charge -q at (0, 0, d/2)
    # Compute multipole at origin
    # We should get [0, 0, 0, q*d, ...] as the multipole moment
    d = 1.0
    q = 5.0
    x1 = np.array([d, 0.0, 0.0])
    x2 = np.array([-d, 0.0, 0.0])
    O = np.array([0.0, 0.0, 0.0])

    for Order in range(1, fmm.FMMGEN_MAXORDER):
        Msize = Nterms(Order)
        S_px = np.zeros(Msize)
        S_px[0] = q
        S_mx = np.zeros(Msize)
        S_mx[0] = -q

        M = np.zeros(Msize)

        fmm.M2M(*(O - x1), S_px, M, Order)
        fmm.M2M(*(O - x2), S_mx, M, Order)

        print(M)

        assert M[0] == 0.0
        assert M[1] == -2 * q * d
        assert M[2] == 0.0
        assert M[3] == 0.0


def test_linear_quadrupole():
    print("test_linear_quadrupole")
    source_order = 0
    order = source_order + TOTALORDER
    cse = True
    atomic = True
    precision = "double"

    fmmgen.generate_code(
        order,
        "LinearQuadrupole",
        precision=precision,
        CSE=cse,
        cython=True,
        potential=True,
        field=True,
        source_order=source_order,
        atomic=atomic,
        minpow=5,
        harmonic_derivs=True,
        language="c",
    )

    import pyximport

    pyximport.install()
    import LinearQuadrupole_wrap as fmm

    # Superimpose two dipoles, s.t. they have
    # opposite orientations, and the +q charges
    # sit on top of each other.

    # Locate charge -q at  (-d, 0, 0)
    # Locate charge -q at  (+d, 0, 0)
    # Locate charge +2q at (0, 0, 0)

    # Compute multipole at origin
    # We should get [0, 0, 0, 0, ...] as the multipole moment
    d = 2.0
    q = 5.0
    x1 = np.array([-d, 0.0, 0.0])
    x2 = np.array([d, 0.0, 0.0])
    O = np.array([0.0, 0.0, 0.0])

    for Order in range(2, fmm.FMMGEN_MAXORDER):
        Msize = Nterms(Order)
        S_mx = np.zeros(Msize)
        S_mx[0] = -q
        S_px = np.zeros(Msize)
        S_px[0] = -q
        S_O = np.zeros(Msize)
        S_O[0] = 2 * q

        M = np.zeros(Msize)

        fmm.M2M(*(O - x1), S_mx, M, Order)
        print(M)
        fmm.M2M(*(O - x2), S_px, M, Order)
        fmm.M2M(*(O - O), S_O, M, Order)
        print(M)

        assert M[0] == 0.0
        assert M[1] == 0.0
        assert M[2] == 0.0
        assert M[3] == 0.0
        assert M[4] == -2 * d * q


def test_quadrupole_two_dipoles():
    print("test_quadrupole_two_dipoles")
    source_order = 0
    order = source_order + TOTALORDER
    cse = True
    atomic = True
    precision = "double"

    fmmgen.generate_code(
        order,
        "QuadrupoleTwoDipoles",
        precision=precision,
        CSE=cse,
        cython=True,
        potential=True,
        field=True,
        source_order=source_order,
        atomic=atomic,
        minpow=5,
        harmonic_derivs=True,
        language="c",
    )

    import pyximport

    pyximport.install()
    import LinearQuadrupole_wrap as fmm

    # Superimpose two dipoles, s.t. they have
    # opposite orientations, and the +q charges
    # sit on top of each other.

    # Locate dipole (mux, 0, 0)  at  (-d/2, 0, 0)
    # Locate charge (-mux, 0, 0) at  (+d/2, 0, 0)

    # Compute multipole at origin
    # We should get [0, 0, 0, 0.0, ...] as the multipole moment
    d = 2.0
    mux = -5.0
    xp = np.array([d, 0.0, 0.0])
    xm = np.array([-d, 0.0, 0.0])
    O = np.array([0.0, 0.0, 0.0])

    for Order in range(2, fmm.FMMGEN_MAXORDER):
        print(f"Quadrupole Test {Order}")
        Msize = Nterms(Order)
        S_px = np.zeros(Msize)
        S_px[1] = -mux
        S_mx = np.zeros(Msize)
        S_mx[1] = mux

        M = np.zeros(Msize)

        fmm.M2M(*(O - xp), S_px, M, Order)
        # print(f'M after dipole ({-mux}, 0, 0) at({}/2, 0, 0) = {M}')
        fmm.M2M(*(O - xm), S_mx, M, Order)
        print(M)
        assert M[0] == 0.0
        assert M[1] == 0.0
        assert M[2] == 0.0
        assert M[3] == 0.0
        assert M[4] == 2 * mux * d


if __name__ == "__main__":
    test_linear_dipole()
    test_linear_quadrupole()
    test_quadrupole_two_dipoles()
