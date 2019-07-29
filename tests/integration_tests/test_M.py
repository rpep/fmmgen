import fmmgen
import numpy as np
from fmmgen.utils import Nterms
import pyximport
pyximport.install()

def test_M2M_operator_0():
    source_order = 0
    order = source_order + 3
    cse = True
    atomic = True
    precision='double'
    
    fmmgen.generate_code(order, "testM2M0",
                         precision=precision,
                         CSE=cse,
                         cython_wrapper=True,
                         potential=True,
                         field=True,
                         source_order=source_order,
                         atomic=atomic, minpow=5,
                         harmonic_derivs=True,
                         language='c')

    import testM2M0_wrap as fmm
    # Locate source S at x1
    # Compute multipoles at z1
    # Compute multipoles using shift at z2 from those at z1
    x1 = np.array([5.0, 5.0, 5.0])
    z1 = np.array([4.0, 4.0, 4.0])
    z2 = np.array([0.0, 0.0, 0.0])

    S = np.zeros(fmm.FMMGEN_SOURCESIZE)
    F = np.zeros(fmm.FMMGEN_OUTPUTSIZE)
    S[:] = 1.0

    for Order in range(fmm.FMMGEN_MINORDER, fmm.FMMGEN_MAXORDER):
        Msize = Nterms(Order)
        S_tmp = np.zeros(Msize)
        S_tmp[:fmm.FMMGEN_SOURCESIZE] = S
        print(S_tmp)

        M_z1 = np.zeros(Msize)
        M_z2_noshift = np.zeros(Msize)
        M_z2_shift = np.zeros(Msize)

        # x1 - z1 for P2M
        fmm.M2M(*(z1 - x1), S_tmp, M_z1, Order)
        fmm.M2M(*(z2 - x1), S_tmp, M_z2_noshift, Order)
        fmm.M2M(*(z2 - z1 ), M_z1, M_z2_shift, Order)

        print(f'Order = {Order}')
        print(f'   S          = {S_tmp}\n')
        print(f'   M_z1       = {M_z1}\n')
        print(f'   z2 noshift = {M_z2_noshift}\n')
        print(f'   z2 shift = {M_z2_shift}\n')
        np.testing.assert_equal(M_z2_shift, M_z2_noshift)

def test_M2M_operator_1():
    source_order = 1
    order = source_order + 3
    cse = True
    atomic = True
    precision='double'

    fmmgen.generate_code(order, "test_M2M1",
                         precision=precision,
                         CSE=cse,
                         cython_wrapper=True,
                         potential=True,
                         field=True,
                         source_order=source_order,
                         atomic=atomic, minpow=5,
                         harmonic_derivs=True,
                         language='c')

    import testM2M1_wrap as fmm
    # Locate source S at x1
    # Compute multipoles at z1
    # Compute multipoles using shift at z2
    x1 = np.array([5.0, 5.0, 5.0])
    z1 = np.array([4.0, 4.0, 4.0])
    z2 = np.array([0.0, 0.0, 0.0])

    S = np.zeros(fmm.FMMGEN_SOURCESIZE)
    F = np.zeros(fmm.FMMGEN_OUTPUTSIZE)
    S[:] = 1.0

    for Order in range(fmm.FMMGEN_MINORDER, fmm.FMMGEN_MAXORDER):
        Msize = Nterms(Order)
        S_tmp = np.zeros(Msize)
        S_tmp[:fmm.FMMGEN_SOURCESIZE] = S

        M_z1 = np.zeros(Msize)
        M_z2_noshift = np.zeros(Msize)
        M_z2_shift = np.zeros(Msize)

        # x1 - z1 for P2M
        fmm.M2M(*(z1 - x1), S_tmp, M_z1, Order)
        fmm.M2M(*(z2 - x1), S_tmp, M_z2_noshift, Order)
        fmm.M2M(*(z2 - z1 ), M_z1, M_z2_shift, Order)
        print(f'Order = {Order}')
        print(f'   S          = {S_tmp}\n')
        print(f'   M_z1       = {M_z1}\n')
        print(f'   z2 noshift = {M_z2_noshift}\n')
        print(f'   z2 shift = {M_z2_shift}\n')
        np.testing.assert_equal(M_z2_shift, M_z2_noshift)
