import fmmgen
import numpy as np
from fmmgen.utils import Nterms
import logging

logging.getLogger("fmmgen")

TOTALORDER = 4


def test_monopole_M2P():
    source_order = 0
    order = source_order + TOTALORDER
    cse = True
    atomic = True
    precision = "double"

    fmmgen.generate_code(
        order,
        "MonopoleOrigin",
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
    import MonopoleOrigin_wrap as fmm

    # Locate charge +q at (0, 0, d/2)
    # Locate charge -q at (0, 0, d/2)
    # Compute multipole at origin
    # We should get [0, 0, 0, q*d, ...] as the multipole moment
    d = 2.0
    q = 5.0
    x1 = np.array([d, 0.0, 0.0])
    x2 = np.array([20.0, 0.0, 0.0])
    dx = x2 - x1
    O = np.array([0.0, 0.0, 0.0])

    Rinv = 1 / np.linalg.norm(dx)
    F_python = np.array(
        [
            q * Rinv,
            q * dx[0] * (Rinv * Rinv * Rinv),
            q * dx[1] * (Rinv * Rinv * Rinv),
            q * dx[2] * (Rinv * Rinv * Rinv),
        ]
    )

    print(f"F_python = {F_python}")

    # Note: order 0 has no generated operator (generate_code(N) emits 1..N-1),
    # so the dispatcher would leave F untouched. Start at order 1.
    F0_last = np.inf
    F1_last = np.inf
    for Order in range(1, fmm.FMMGEN_MAXORDER):
        Msize = Nterms(Order)
        S = np.zeros(Msize)
        S[0] = q
        M = np.zeros(Msize)
        F = np.zeros(4)
        fmm.M2M(*(O - x1), S, M, Order)
        print(M)
        fmm.M2P(*(x2 - O), M, F, Order)
        print(f"Order = {Order}, F = {F}")

        F0_err = np.abs((F[0] - F_python[0]) / F_python[0])
        F1_err = np.abs((F[1] - F_python[1]) / F_python[1])

        # The error must strictly decrease with expansion order. The previous
        # version of this test compared against a constant 100 that was never
        # updated, so a completely wrong (or zero) field passed.
        assert F0_err < F0_last, f"potential not converging at order {Order}: {F0_err} !< {F0_last}"
        assert F1_err < F1_last, f"field not converging at order {Order}: {F1_err} !< {F1_last}"
        assert F[2] == 0
        assert F[3] == 0

        F0_last = F0_err
        F1_last = F1_err

    # Absolute accuracy at the highest available order.
    assert F0_last < 1e-3, f"potential error too large at highest order: {F0_last}"
    assert F1_last < 5e-3, f"field error too large at highest order: {F1_last}"


def test_monopole_M2P_offaxis():
    """M2P field for a geometry with all three components nonzero.

    test_monopole_M2P is collinear, so Fy and Fz are identically zero there and
    a sign or indexing error in those components cannot be detected.
    """
    import pyximport

    pyximport.install()
    import MonopoleOrigin_wrap as fmm

    q = 3.0
    xs = np.array([0.4, -0.3, 0.5])
    xt = np.array([7.0, 5.0, -4.0])
    O = np.zeros(3)

    d = xt - xs
    rn = np.linalg.norm(d)
    ref = np.array([q / rn, *(q * d / rn**3)])

    last = np.inf
    for Order in range(1, fmm.FMMGEN_MAXORDER):
        Msize = Nterms(Order)
        S = np.zeros(Msize)
        S[0] = q
        M = np.zeros(Msize)
        F = np.zeros(4)
        fmm.M2M(*(O - xs), S, M, Order)
        fmm.M2P(*(xt - O), M, F, Order)
        err = np.linalg.norm(F[1:] - ref[1:]) / np.linalg.norm(ref[1:])
        print(f"Order = {Order}, field rel err = {err:.3e}")
        assert err < last, f"field not converging at order {Order}: {err} !< {last}"
        last = err

    assert last < 1e-3, f"field error too large at highest order: {last}"


def test_P2P_field_monopole():
    """P2P is exact, so both potential and field must match analytically.

    This is a regression test: the P2P field operators were previously generated
    as identically zero.
    """
    import pyximport

    pyximport.install()
    import MonopoleOrigin_wrap as fmm

    q = 2.5
    d = np.array([1.5, -2.0, 3.5])
    rn = np.linalg.norm(d)
    ref = np.array([q / rn, *(q * d / rn**3)])

    S = np.zeros(fmm.FMMGEN_SOURCESIZE)
    S[0] = q
    F = np.zeros(4)
    fmm.P2P(*d, S, F)

    print(f"P2P F = {F}, reference = {ref}")
    assert np.linalg.norm(F[1:]) > 0, "P2P field is identically zero"
    np.testing.assert_allclose(F, ref, rtol=1e-12)
