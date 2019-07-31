import fmmgen
import numpy as np
from fmmgen.utils import Nterms
import pytest
import pyximport
pyximport.install()


def test_P2P_operator_0():
    source_order = 0
    order = source_order + 1
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
    x2 = np.array([4.0, 4.0, 4.0])

    S = np.zeros(fmm.FMMGEN_SOURCESIZE)
    F_python = np.zeros(fmm.FMMGEN_OUTPUTSIZE)
    F_c = np.zeros(fmm.FMMGEN_OUTPUTSIZE)

    S[:] = 1.0

    dx = x2 - x1
    R = np.linalg.norm(dx)
    
    F_python = np.array([S[0] / R, S[0]/R**3, S[0]/R**3, S[0]/R**3])
    
    fmm.P2P(*dx, S, F_c)
    np.testing.assert_equal(F_python, F_c)



def test_P2P_operator_1():
    source_order = 1
    order = source_order + 1
    cse = True
    atomic = True
    precision='double'

    fmmgen.generate_code(order, "testM2M1",
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
    # Compute multipoles using shift at z2 from those at z1
    x1 = np.array([5.0, 5.0, 5.0])
    x2 = np.array([4.0, 4.0, 4.0])

    S = np.zeros(fmm.FMMGEN_SOURCESIZE)
    F_python = np.zeros(fmm.FMMGEN_OUTPUTSIZE)
    F_c = np.zeros(fmm.FMMGEN_OUTPUTSIZE)

    S[:] = 1.0
    print(S)
    dx = x2 - x1
    R = np.linalg.norm(dx)
    fmm.P2P(*dx, S, F_c)
    
    dx /= R
    F_python = np.array([S.dot(dx)/R**2, 0, 0, 0])
    F_python[1:] = (3 * S.dot(dx) * dx - S) / R**3 
    

    np.testing.assert_equal(F_python, F_c)
