import fmmgen
import numpy as np
from fmmgen.utils import Nterms
import logging
logging.getLogger('fmmgen')

TOTALORDER = 4

def test_monopole_M2P():
    source_order = 0
    order = source_order + TOTALORDER
    cse = True
    atomic = True
    precision='double'

    fmmgen.generate_code(order, "MonopoleOrigin",
                         precision=precision,
                         CSE=cse,
                         cython=True,
                         potential=True,
                         field=True,
                         source_order=source_order,
                         atomic=atomic, minpow=5,
                         harmonic_derivs=True,
                         language='c')

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
    F_python = np.array([q*Rinv,
                         q*dx[0]*(Rinv*Rinv*Rinv),
                         q*dx[1]*(Rinv*Rinv*Rinv),
                         q*dx[2]*(Rinv*Rinv*Rinv)])

    print(f'F_python = {F_python}')

    F0_last = 100
    F1_last = 100
    for Order in range(0, fmm.FMMGEN_MAXORDER):
        Msize = Nterms(Order)
        S = np.zeros(Msize)
        S[0] = q
        M = np.zeros(Msize)
        F = np.zeros(4)
        fmm.M2M(*(O-x1), S, M, Order)
        print(M)
        fmm.M2P(*(x2-O), M, F, Order)
        print(f'Order = {Order}, F = {F}')

        assert np.abs((F[0] - F_python[0])/F_python[0]) < F0_last
        assert np.abs((F[1] - F_python[1])/F_python[1]) < F1_last
        assert F[2] == 0
        assert F[3] == 0
