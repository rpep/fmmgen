"""
Generated kernel-generic operators, compiled with gcc and driven through
ctypes: M2L/L2P against direct summation with geometric convergence in p, and
the Kernel(1/r) path against the original Laplace operators.

Needs gcc (and g++ for the C++ smoke test); skipped without.
"""
import os

import numpy as np
import pytest
import sympy as sp

from fmmgen import generate_code
from fmmgen.kernel import multi_indices, nterms
from fmmgen.utils import generate_mappings

from kernel_util import (HAVE_GCC, HAVE_GXX, KERNELS, Lib, far_field, rel_err)

pytestmark = pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")

ORDERS = [1, 3, 5]


@pytest.fixture(scope="module", params=sorted(KERNELS))
def compiled(request, tmp_path_factory):
    k = KERNELS[request.param]()
    return request.param, Lib(k, tmp_path_factory.mktemp(request.param), ORDERS)


def test_far_field_converges(compiled):
    """P2M -> M2L -> L2P vs direct sum; error falls geometrically in p."""
    name, lib = compiled
    errs = {}
    for p in ORDERS:
        rng = np.random.default_rng(11)
        # the Gaussian has width ~1, so its clusters must be small against that
        geo = dict(dist=1.5, radius=0.15) if name == "gaussian" else {}
        phi, F, pr, Fr = far_field(lib, p, rng, "A", **geo)
        errs[p] = (rel_err(phi, pr), rel_err(F, Fr))
    for i in (0, 1):
        # strictly improving, and by a clear factor from p=1 to p=5
        assert errs[3][i] < errs[1][i]
        assert errs[5][i] < errs[3][i] or errs[5][i] < 1e-9
        assert errs[5][i] < 0.05 * errs[1][i], (name, errs)
    assert errs[5][0] < 5e-3, (name, errs)


def test_translations_exact(compiled):
    """M2M and L2L are exact on the truncated expansions: the full route
    P2M(children) -> M2M -> M2L -> L2L -> L2P must equal P2M(parent) -> M2L
    -> L2P to rounding, on identical points."""
    name, lib = compiled
    d = lib.dim
    p = 3
    N = nterms(d, p)
    rng = np.random.default_rng(3)
    cs, ct = np.zeros(d), np.full(d, 3.0 / np.sqrt(d))
    c1, c2 = cs + 0.2, cs - 0.15
    srcs = [c1 + 0.1 * rng.standard_normal(d), c2 + 0.1 * rng.standard_normal(d)]
    q = [0.7, -1.3]
    # route 1: direct at the parent
    M = np.zeros(N)
    for s, qq in zip(srcs, q):
        lib.call(f"S2M_{p}", cs - s, np.array([qq]), M)
    # route 2: children then M2M
    M2 = np.zeros(N)
    for c, s, qq in zip((c1, c2), srcs, q):
        Mk = np.zeros(N)
        lib.call(f"S2M_{p}", c - s, np.array([qq]), Mk)
        lib.call(f"M2M_{p}", cs - c, Mk, M2)
    assert np.allclose(M, M2, rtol=1e-11, atol=1e-12)
    L = np.zeros(N)
    lib.call(f"M2L_{p}", ct - cs, M, L)
    child = ct + 0.1
    Lc = np.zeros(N)
    lib.call(f"L2L_{p}", child - ct, L, Lc)
    pt = child + 0.05 * rng.standard_normal(d)
    F1, F2 = np.zeros(1 + d), np.zeros(1 + d)
    lib.call(f"L2P_{p}", pt - ct, L, F1)
    lib.call(f"L2P_{p}", pt - child, Lc, F2)
    assert np.allclose(F1, F2, rtol=1e-10, atol=1e-12)


def test_p2p_and_m2p(compiled):
    name, lib = compiled
    k, d = lib.kernel, lib.dim
    rng = np.random.default_rng(2)
    s = rng.standard_normal(d) * 0.1
    t = s + np.array([1.1, -0.7, 0.9][:d])
    G, grads = k.lambdify()
    F = np.zeros(1 + d)
    lib.call("P2P", t - s, np.array([2.0]), F)
    assert np.isclose(F[0], 2.0 * G(*(t - s)), rtol=1e-12)
    for j in range(d):
        assert np.isclose(F[1 + j], -2.0 * grads[j](*(t - s)), rtol=1e-10, atol=1e-13)
    # M2P with a multipole that is exactly a point charge at the centre
    p = 5
    M = np.zeros(nterms(d, p))
    lib.call(f"S2M_{p}", np.zeros(d), np.array([2.0]), M)
    F = np.zeros(1 + d)
    lib.call(f"M2P_{p}", t - 0 * s, M, F)
    ref = np.zeros(1 + d)
    ref[0] = 2.0 * G(*t)
    ref[1:] = [-2.0 * g(*t) for g in grads]
    assert np.allclose(F, ref, rtol=1e-9, atol=1e-12)


def test_matches_original_laplace(tmp_path):
    """Kernel(1/r) reproduces the original generator's operators."""
    p = 4
    old = tmp_path / "old"
    new = tmp_path / "new"
    old.mkdir()
    new.mkdir()
    cwd = os.getcwd()
    os.chdir(old)
    try:
        generate_code(p + 1, "oldops", CSE=True, include_dir=str(old), src_dir=str(old), language="c")
    finally:
        os.chdir(cwd)
    import subprocess

    subprocess.run(["gcc", "-O1", "-shared", "-fPIC", "-o", str(old / "o.so"), str(old / "oldops.c"), "-lm"],
                   check=True)
    import ctypes

    olib = ctypes.CDLL(str(old / "o.so"))
    nlib = Lib(KERNELS["laplace"](), new, [p])

    class Old:
        dim = 3

        def call(self, fname, coords, *arrays):
            dp = ctypes.POINTER(ctypes.c_double)
            f = getattr(olib, fname)
            f.argtypes = [ctypes.c_double] * 3 + [dp] * len(arrays)
            f.restype = None
            f(*[float(c) for c in coords], *[a.ctypes.data_as(dp) for a in arrays])

    O = Old()
    rng = np.random.default_rng(8)
    N = nterms(3, p)
    xv = rng.standard_normal(3)
    for fn, nin, nout in ((f"M2M_{p}", N, N), (f"L2L_{p}", N, N), (f"M2L_{p}", N, N),
                          (f"L2P_{p}", N, 4), (f"M2P_{p}", N, 4)):
        v = rng.standard_normal(nin)
        xx = xv * (3.0 if "M2L" in fn or "M2P" in fn else 0.4)
        a, b = np.zeros(nout), np.zeros(nout)
        O.call(fn, xx, v, a)
        nlib.call(fn, xx, v, b)
        assert np.allclose(a, b, rtol=1e-10, atol=1e-11), fn
    a, b = np.zeros(N), np.zeros(N)
    O.call(f"S2M_{p}", xv, np.array([1.7]), a)
    nlib.call(f"S2M_{p}", xv, np.array([1.7]), b)
    assert np.allclose(a, b, rtol=1e-12)
    a, b = np.zeros(4), np.zeros(4)
    O.call("P2P", xv, np.array([1.7]), a)
    nlib.call("P2P", xv, np.array([1.7]), b)
    assert np.allclose(a, b, rtol=1e-10)


def test_index_ordering_matches_legacy():
    sym = sp.symbols("x y z")
    for p in (1, 3, 5):
        legacy, _ = generate_mappings(p, sym, "grevlex")
        assert legacy == multi_indices(3, p)


@pytest.mark.skipif(not HAVE_GXX, reason="g++ not available")
def test_cxx_output_compiles(tmp_path):
    """The same generator through the C++ printer."""
    Lib(KERNELS["yukawa"](), tmp_path, [2], language="c++")


def test_potential_only_and_field_only(tmp_path):
    k = KERNELS["plummer"]()
    lib = Lib(k, tmp_path, [3], potential=True, field=False)
    rng = np.random.default_rng(1)
    M = np.zeros(nterms(3, 3))
    lib.call("S2M_3", rng.standard_normal(3), np.array([1.0]), M)
    F = np.zeros(1)
    lib.call("P2P", np.array([1.0, 2.0, 0.5]), np.array([1.0]), F)
    assert np.isclose(F[0], 1 / np.sqrt(1 + 4 + 0.25 + 0.09))
