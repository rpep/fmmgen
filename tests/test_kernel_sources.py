"""Multipole sources (dipoles, quadrupoles, octupoles) on the kernel path:
direct summation, M2M/L2L exactness and the legacy source_order operators."""
import ctypes
import os
import subprocess

import numpy as np
import pytest

from fmmgen import generate_code
from fmmgen.kernel import nterms, shell_indices, source_indices
from kernel_util import HAVE_GCC, KERNELS, Lib, direct_s, far_field_s, rel_err

pytestmark = pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")

CASES = [("laplace", 1), ("laplace", 2), ("laplace", 3), ("yukawa", 1), ("yukawa", 2),
         ("plummer", 1), ("plummer", 2), ("aniso", 1), ("log2d", 2)]
ORDERS = {1: [2, 4, 6], 2: [3, 5, 7], 3: [4, 6, 8]}


@pytest.fixture(scope="module", params=CASES, ids=lambda c: f"{c[0]}-s{c[1]}")
def compiled(request, tmp_path_factory):
    name, s = request.param
    k = KERNELS[name]()
    return name, s, Lib(k, tmp_path_factory.mktemp(f"{name}{s}"), ORDERS[s], source_order=s)


def test_p2p_direct(compiled):
    name, s, lib = compiled
    d = lib.dim
    sh = shell_indices(d, s)
    rng = np.random.default_rng(4)
    y_ = rng.standard_normal(d) * 0.1
    t = y_ + np.array([1.1, -0.7, 0.9][:d])
    S = rng.standard_normal((1, len(sh)))
    F = np.zeros(1 + d)
    lib.call("P2P", t - y_, np.ascontiguousarray(S[0]), F)
    phi, Fr = direct_s(lib.kernel, t[None, :], y_[None, :], S, s)
    assert np.isclose(F[0], phi[0], rtol=1e-10, atol=1e-12)
    assert np.allclose(F[1:], Fr[0], rtol=1e-9, atol=1e-11)


def test_far_field_converges(compiled):
    name, s, lib = compiled
    errs = []
    for p in ORDERS[s]:
        phi, F, pr, Fr = far_field_s(lib, p, s, np.random.default_rng(7))
        errs.append((rel_err(phi, pr), rel_err(F, Fr)))
    assert errs[1][0] < errs[0][0] and errs[2][0] < errs[1][0] * 1.5
    assert errs[2][0] < 0.05 * errs[0][0], (name, s, errs)
    assert errs[2][0] < 2e-2, (name, s, errs)


def test_translations_exact(compiled):
    """children -> M2M -> M2L -> L2L -> L2P equals the parent route."""
    name, s, lib = compiled
    d, p = lib.dim, ORDERS[s][0]
    NM, NL = len(source_indices(d, p, s)), nterms(d, p - s)
    sh = shell_indices(d, s)
    rng = np.random.default_rng(3)
    cs, ct = np.zeros(d), np.full(d, 3.0 / np.sqrt(d))
    c1, c2 = cs + 0.2, cs - 0.15
    srcs = [c1 + 0.1 * rng.standard_normal(d), c2 + 0.1 * rng.standard_normal(d)]
    S = rng.standard_normal((2, len(sh)))
    M, M2 = np.zeros(NM), np.zeros(NM)
    for c, y_, Sj in zip((c1, c2), srcs, S):
        lib.call(f"S2M_{p}", cs - y_, np.ascontiguousarray(Sj), M)
        Mk = np.zeros(NM)
        lib.call(f"S2M_{p}", c - y_, np.ascontiguousarray(Sj), Mk)
        lib.call(f"M2M_{p}", cs - c, Mk, M2)
    assert np.allclose(M, M2, rtol=1e-10, atol=1e-12)
    L = np.zeros(NL)
    lib.call(f"M2L_{p}", ct - cs, M, L)
    child = ct + 0.1
    Lc = np.zeros(NL)
    lib.call(f"L2L_{p}", child - ct, L, Lc)
    pt = child + 0.05 * rng.standard_normal(d)
    F1, F2 = np.zeros(1 + d), np.zeros(1 + d)
    lib.call(f"L2P_{p}", pt - ct, L, F1)
    lib.call(f"L2P_{p}", pt - child, Lc, F2)
    assert np.allclose(F1, F2, rtol=1e-9, atol=1e-11)


def test_m2p_matches_direct_expansion(compiled):
    name, s, lib = compiled
    d, p = lib.dim, ORDERS[s][-1]
    sh = shell_indices(d, s)
    rng = np.random.default_rng(9)
    y_ = rng.standard_normal(d) * 0.05
    S = rng.standard_normal((1, len(sh)))
    M = np.zeros(len(source_indices(d, p, s)))
    lib.call(f"S2M_{p}", -y_, np.ascontiguousarray(S[0]), M)
    t = np.array([2.0, -1.3, 1.1][:d])
    F = np.zeros(1 + d)
    lib.call(f"M2P_{p}", t, M, F)
    phi, Fr = direct_s(lib.kernel, t[None, :], y_[None, :], S, s)
    assert np.isclose(F[0], phi[0], rtol=2e-3)
    assert np.allclose(F[1:], Fr[0], rtol=2e-2, atol=1e-4)


@pytest.mark.parametrize("s", [1, 2, 3])
def test_matches_legacy_source_order(s, tmp_path):
    """Kernel(1/r) with source_order=s reproduces the legacy operators."""
    p = s + 3
    old, new = tmp_path / "old", tmp_path / "new"
    old.mkdir()
    new.mkdir()
    cwd = os.getcwd()
    os.chdir(old)
    try:
        generate_code(p + 1, "oldops", CSE=True, source_order=s, include_dir=str(old), src_dir=str(old),
                      language="c")
    finally:
        os.chdir(cwd)
    subprocess.run(["gcc", "-O1", "-shared", "-fPIC", "-o", str(old / "o.so"), str(old / "oldops.c"), "-lm"],
                   check=True)
    olib = ctypes.CDLL(str(old / "o.so"))
    nlib = Lib(KERNELS["laplace"](), new, [p], source_order=s)
    dp = ctypes.POINTER(ctypes.c_double)

    def ocall(fname, coords, *arrays):
        f = getattr(olib, fname)
        f.argtypes = [ctypes.c_double] * 3 + [dp] * len(arrays)
        f.restype = None
        f(*[float(c) for c in coords], *[a.ctypes.data_as(dp) for a in arrays])

    rng = np.random.default_rng(8)
    NM, NL, sh = len(source_indices(3, p, s)), nterms(3, p - s), len(shell_indices(3, s))
    xv = rng.standard_normal(3)
    Sv = rng.standard_normal(sh)

    def both(fname, xx, nout, *ins):
        a, b = np.zeros(nout), np.zeros(nout)
        ocall(fname, xx, *ins, a)
        nlib.call(fname, xx, *ins, b)
        assert np.allclose(a, b, rtol=1e-9, atol=1e-10), (fname, s, np.max(np.abs(a - b)))

    both(f"S2M_{p}", xv * 0.4, NM, Sv)
    both(f"M2M_{p}", xv * 0.4, NM, rng.standard_normal(NM))
    both(f"M2L_{p}", xv * 3.0, NL, rng.standard_normal(NM))
    both(f"L2L_{p}", xv * 0.4, NL, rng.standard_normal(NL))
    both(f"L2P_{p}", xv * 0.4, 4, rng.standard_normal(NL))
    both(f"M2P_{p}", xv * 3.0, 4, rng.standard_normal(NM))
    both("P2P", xv, 4, Sv)
