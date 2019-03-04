# Generate the code
import example 
# Import with Cython
import pyximport
pyximport.install()

import numpy as np
import fmm_wrap
from fmmgen.utils import Nterms

def test_fmm_wrap_P2M0():
    M = np.zeros(Nterms(0))
    fmm_wrap.P2M0(0, 0, 0, 1.0, M)
    assert M[0] == 1.0


def test_fmm_wrap_P2M1():
    M = np.zeros(Nterms(1))
    fmm_wrap.P2M1(2.0, 2.0, 2.0, 1.0, M)
    assert M[0] == 1.0
    assert M[1] == -2.0
    assert M[2] == -2.0
    assert M[3] == -2.0

def test_fmm_wrap_P2M2():
    M = np.zeros(Nterms(2))
    fmm_wrap.P2M2(2.0, 2.0, 2.0, 1.0, M)
    assert M[0] == 1.0   # q
    assert M[1] == -2.0  # -q*x
    assert M[2] == -2.0  # -q*y
    assert M[3] == -2.0  # -q*z
    assert M[4] == 2.0   # 0.5 * q*x**2
    assert M[5] == 4.0   # q*x*y
    assert M[6] == 4.0   # q*x*z
    assert M[7] == 2.0   # 0.5 * q*y**2
    assert M[8] == 4.0   # q*y*z
    assert M[9] == 2.0   # 0.5 * q*z**2
