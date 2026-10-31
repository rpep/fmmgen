"""Spherical-harmonic Laplace operators: numpy reference (fmmgen.sh_math) and the
generated C code (fmmgen.sh), checked against direct summation and against each
other."""
import shutil

import numpy as np
import pytest

from fmmgen import sh_math as sm
from fmmgen.sh import generate_sh_code, rot_size

from sh_util import SHLib, cluster, direct, far_field_chain, rel_err

pytestmark = pytest.mark.skipif(shutil.which("gcc") is None, reason="needs gcc")

CS = np.zeros(3)
CT = np.array([3.0, 1.0, -2.0])


def _problem(rng, nsrc=20, ntgt=5, half=0.5):
    src = cluster(rng, CS, nsrc, half)
    q = rng.standard_normal(nsrc)
    tgt = cluster(rng, CT, ntgt, half)
    return src, q, tgt


def _multipole(p, src, q, centre=CS):
    M = np.zeros(sm.sizes(p)[0])
    for s, qq in zip(src, q):
        sm.S2M(p, *(centre - s), np.array([qq]), M)
    return M


def _err(p, L, tgt, ct, src, q):
    F = np.zeros((len(tgt), 4))
    for i, t in enumerate(tgt):
        sm.L2P(p, *(t - ct), L, F[i])
    phi, fld = direct(tgt, src, q)
    return rel_err(F[:, 0], phi), rel_err(F[:, 1:], fld)


def test_sizes():
    for p in range(0, 17):
        assert sm.sizes(p) == ((p + 1) ** 2, (p + 1) ** 2)


def test_solid_harmonic_expansion_of_one_over_r():
    """1/|r - r'| = sum conj(R_n^m(r')) I_n^m(r) in the flat real layout."""
    p = 14
    rng = np.random.default_rng(1)
    rp = rng.standard_normal(3) * 0.3
    r = np.array([2.0, 1.0, -1.5])
    R, I = sm.solid_R(p, *rp), sm.solid_I(p, *r)
    s = R[0] * I[0]
    for n in range(1, p + 1):
        s += R[sm.ire(n, 0)] * I[sm.ire(n, 0)]
        for m in range(1, n + 1):
            s += 2 * (R[sm.ire(n, m)] * I[sm.ire(n, m)] + R[sm.iim(n, m)] * I[sm.iim(n, m)])
    assert abs(s - 1 / np.linalg.norm(r - rp)) < 1e-10


def test_p2p_matches_direct():
    rng = np.random.default_rng(2)
    src = rng.standard_normal((7, 3))
    q = rng.standard_normal(7)
    t = np.array([2.0, -1.0, 0.5])
    F = np.zeros(4)
    for s, qq in zip(src, q):
        sm.P2P(*(t - s), np.array([qq]), F)
    phi, fld = direct(t[None], src, q)
    assert abs(F[0] - phi[0]) < 1e-13 and np.allclose(F[1:], fld[0], atol=1e-13)


@pytest.mark.parametrize("variant", ["plain", "rotation"])
def test_far_field_converges(variant):
    """Potential and field error against direct summation shrink geometrically in p."""
    rng = np.random.default_rng(3)
    src, q, tgt = _problem(rng)
    errs = []
    for p in (2, 4, 8, 12, 16):
        M = _multipole(p, src, q)
        L = np.zeros_like(M)
        sm.M2L(p, *(CT - CS), M, L, variant)
        errs.append(_err(p, L, tgt, CT, src, q))
    phi_e = [e[0] for e in errs]
    fld_e = [e[1] for e in errs]
    assert all(b < a for a, b in zip(phi_e, phi_e[1:]))
    assert all(b < a for a, b in zip(fld_e, fld_e[1:]))
    assert phi_e[-1] < 1e-8 and fld_e[-1] < 1e-6


def test_m2p_converges():
    rng = np.random.default_rng(4)
    src, q, tgt = _problem(rng)
    prev = 1.0
    for p in (2, 6, 10, 14):
        M = _multipole(p, src, q)
        F = np.zeros((len(tgt), 4))
        for i, t in enumerate(tgt):
            sm.M2P(p, *(t - CS), M, F[i])
        phi, fld = direct(tgt, src, q)
        e = max(rel_err(F[:, 0], phi), rel_err(F[:, 1:], fld))
        assert e < prev
        prev = e
    assert prev < 1e-6


@pytest.mark.parametrize("p", [3, 8, 14])
def test_m2m_is_exact(p):
    """S2M about the parent centre equals M2M of S2M about the child centre (degree p, exactly)."""
    rng = np.random.default_rng(5)
    src = rng.uniform(-0.4, 0.4, (6, 3))
    q = rng.standard_normal(6)
    child = np.array([0.3, -0.2, 0.25])
    parent = np.array([0.1, 0.0, -0.1])
    Mc = _multipole(p, src, q, child)
    Mp = np.zeros_like(Mc)
    sm.M2M(p, *(parent - child), Mc, Mp)
    assert rel_err(Mp, _multipole(p, src, q, parent)) < 1e-12


@pytest.mark.parametrize("p", [3, 8, 14])
def test_l2l_is_exact(p):
    """L2P about the child centre of L2L(L) equals L2P about the parent centre of L at any point."""
    rng = np.random.default_rng(6)
    L = rng.standard_normal(sm.sizes(p)[0])
    child = np.array([0.2, 0.1, -0.2])
    Ls = np.zeros_like(L)
    sm.L2L(p, *child, L, Ls)
    for _ in range(4):
        t = child + rng.uniform(-0.4, 0.4, 3)
        a, b = np.zeros(4), np.zeros(4)
        sm.L2P(p, *t, L, a)
        sm.L2P(p, *(t - child), Ls, b)
        assert np.allclose(a, b, rtol=1e-11, atol=1e-11 * np.abs(a).max())


def test_field_is_minus_gradient_of_potential():
    rng = np.random.default_rng(7)
    p = 9
    L = rng.standard_normal(sm.sizes(p)[0])
    t = np.array([0.3, -0.2, 0.4])
    F = np.zeros(4)
    sm.L2P(p, *t, L, F)
    h = 1e-6
    for i in range(3):
        e = np.zeros(3)
        e[i] = h
        a, b = np.zeros(4), np.zeros(4)
        sm.L2P(p, *(t + e), L, a)
        sm.L2P(p, *(t - e), L, b)
        assert abs(F[1 + i] + (a[0] - b[0]) / (2 * h)) < 1e-6 * np.abs(F).max()
    M = rng.standard_normal(sm.sizes(p)[0])
    t = np.array([2.0, -1.5, 1.0])
    F = np.zeros(4)
    sm.M2P(p, *t, M, F)
    for i in range(3):
        e = np.zeros(3)
        e[i] = h
        a, b = np.zeros(4), np.zeros(4)
        sm.M2P(p, *(t + e), M, a)
        sm.M2P(p, *(t - e), M, b)
        assert abs(F[1 + i] + (a[0] - b[0]) / (2 * h)) < 1e-6 * np.abs(F).max()


@pytest.mark.parametrize("trunc", ["box", "total"])
def test_plain_and_rotation_m2l_agree(trunc):
    """The rotation-based M2L reproduces the plain double sum to rounding, p = 2..16,
    for random multipoles and directions including both poles."""
    rng = np.random.default_rng(8)
    dirs = [rng.uniform(-3, 3, 3) for _ in range(2)] + [np.array([0, 0, 2.5]), np.array([0, 0, -2.5]),
                                                          np.array([2.0, 0, 0])]
    for p in range(2, 17):
        M = rng.standard_normal(sm.sizes(p)[0])
        for d in dirs:
            a, b = np.zeros_like(M), np.zeros_like(M)
            sm.m2l_plain(p, *d, M, a, trunc)
            sm.m2l_rotation(p, *d, M, b, trunc)
            assert rel_err(a, b) < 1e-11, (p, d)


def test_rotation_matrices_are_orthogonal_and_map_d_to_z():
    rng = np.random.default_rng(9)
    d = rng.standard_normal(3)
    Q, rho = sm.rotation_to_z(*d)
    assert np.allclose(Q @ Q.T, np.eye(3)) and np.allclose(Q @ d, [0, 0, rho])
    for D in sm.ir_rotations(Q, 12):
        assert np.allclose(D @ D.T, np.eye(len(D)), atol=1e-12)


# ---------------------------------------------------------------------------
# Generated C
# ---------------------------------------------------------------------------
@pytest.fixture(scope="module")
def loop_lib(tmp_path_factory):
    return SHLib(tmp_path_factory.mktemp("shloop"), 17, name="shloop", unroll=False)


@pytest.fixture(scope="module")
def unrolled_lib(tmp_path_factory):
    return SHLib(tmp_path_factory.mktemp("shunr"), 6, name="shunr", unroll=True)


def _ops_outputs(lib, p, rng, trunc_total=False):
    N = sm.sizes(p)[0]
    v = rng.uniform(-1, 1, 3) * 2 + np.array([0.3, 0.2, 2.5])
    M = rng.standard_normal(N)
    L0 = rng.standard_normal(N)
    out = {}
    a = np.zeros(N)
    lib.call(f"S2M_{p}", v, np.array([0.7]), a)
    out["S2M"] = a
    for name in ("M2M", "L2L", "M2Lplain", "M2Lrot", "M2L"):
        a = L0.copy()
        lib.call(f"{name}_{p}", v, M, a)
        out[name] = a
    for name in ("L2P", "M2P"):
        a = np.zeros(4)
        lib.call(f"{name}_{p}", v, M, a)
        out[name] = a
    return v, M, L0, out


@pytest.mark.parametrize("p", [1, 2, 5, 9, 16])
def test_loop_c_matches_numpy(loop_lib, p):
    rng = np.random.default_rng(10 + p)
    v, M, L0, out = _ops_outputs(loop_lib, p, rng)
    N = len(M)
    ref = {}
    ref["S2M"] = np.zeros(N)
    sm.S2M(p, *v, np.array([0.7]), ref["S2M"])
    for name, f in (("M2M", sm.M2M), ("L2L", sm.L2L)):
        ref[name] = L0.copy()
        f(p, *v, M, ref[name])
    for name in ("M2Lplain", "M2Lrot", "M2L"):
        ref[name] = L0.copy()
        (sm.m2l_plain if name == "M2Lplain" else sm.m2l_rotation)(p, *v, M, ref[name])
    for name, f in (("L2P", sm.L2P), ("M2P", sm.M2P)):
        ref[name] = np.zeros(4)
        f(p, *v, M, ref[name])
    for k in ref:
        assert rel_err(out[k], ref[k]) < 1e-11, k


def test_loop_c_far_field_chain(loop_lib):
    """S2M -> M2L -> L2P (and the M2M/L2L route) against direct summation, both M2L variants."""
    rng = np.random.default_rng(11)
    for p, tol in ((4, 5e-3), (8, 1e-4), (12, 1e-6), (16, 1e-8)):
        for m2l in ("M2Lplain", "M2Lrot"):
            for route in ("A", "B"):
                phi, F, phi_ref, F_ref = far_field_chain(loop_lib, p, rng, m2l=m2l, route=route)
                assert rel_err(phi, phi_ref) < tol, (p, m2l, route)
                assert rel_err(F, F_ref) < 30 * tol, (p, m2l, route)


def test_loop_c_rotation_setup_apply(loop_lib):
    """Cached-rotation path (setup once, apply many) equals the on-the-fly M2L."""
    rng = np.random.default_rng(12)
    p = 10
    d = np.array([2.0, -3.0, 2.0])
    Dt = np.zeros(rot_size(p))
    loop_lib.call(f"M2Lrot_setup_{p}", d, Dt)
    assert abs(Dt[0] - 1 / np.linalg.norm(d)) < 1e-15
    for _ in range(3):
        M = rng.standard_normal(sm.sizes(p)[0])
        a, b = np.zeros_like(M), np.zeros_like(M)
        loop_lib.apply(f"M2Lrot_apply_{p}", Dt, M, a)
        loop_lib.call(f"M2Lrot_{p}", d, M, b)
        assert rel_err(a, b) < 1e-14


def test_p2p_c(loop_lib):
    rng = np.random.default_rng(13)
    src = rng.standard_normal((9, 3))
    q = rng.standard_normal(9)
    t = np.array([1.5, 2.0, -1.0])
    F = np.zeros(4)
    for s, qq in zip(src, q):
        loop_lib.call("P2P", t - s, np.array([qq]), F)
    phi, fld = direct(t[None], src, q)
    assert abs(F[0] - phi[0]) < 1e-13 and np.allclose(F[1:], fld[0], atol=1e-13)


@pytest.mark.parametrize("p", [1, 3, 5])
def test_unrolled_matches_loop(loop_lib, unrolled_lib, p):
    rng = np.random.default_rng(20 + p)
    v, M, L0, a = _ops_outputs(unrolled_lib, p, rng)
    rng = np.random.default_rng(20 + p)
    _, _, _, b = _ops_outputs(loop_lib, p, rng)
    for k in a:
        assert rel_err(a[k], b[k]) < 1e-12, k


def test_unrolled_rotation_apply_matches_loop(loop_lib, unrolled_lib):
    rng = np.random.default_rng(30)
    p = 4
    d = np.array([-1.0, 2.5, 2.0])
    Dt = np.zeros(rot_size(p))
    unrolled_lib.call(f"M2Lrot_setup_{p}", d, Dt)
    M = rng.standard_normal(sm.sizes(p)[0])
    a, b = np.zeros_like(M), np.zeros_like(M)
    unrolled_lib.apply(f"M2Lrot_apply_{p}", Dt, M, a)
    Dt2 = np.zeros(rot_size(p))
    loop_lib.call(f"M2Lrot_setup_{p}", d, Dt2)
    loop_lib.apply(f"M2Lrot_apply_{p}", Dt2, M, b)
    assert rel_err(a, b) < 1e-13


def test_total_truncation_c(tmp_path):
    lib = SHLib(tmp_path, 7, name="shtot", trunc="total")
    rng = np.random.default_rng(31)
    p = 6
    N = sm.sizes(p)[0]
    v = np.array([1.5, -2.0, 2.5])
    M = rng.standard_normal(N)
    for name, f in (("M2Lplain", sm.m2l_plain), ("M2Lrot", sm.m2l_rotation)):
        a, b = np.zeros(N), np.zeros(N)
        lib.call(f"{name}_{p}", v, M, a)
        f(p, *v, M, b, "total")
        assert rel_err(a, b) < 1e-11


def test_unroll_limit_is_enforced(tmp_path):
    with pytest.raises(ValueError):
        generate_sh_code(12, "x", unroll=True, include_dir=str(tmp_path))


@pytest.mark.skipif(shutil.which("g++") is None, reason="needs g++")
def test_cpp_output_compiles(tmp_path):
    """language='c++' output compiles with g++ -Wall and gives the same M2L as the numpy reference
    (C++ names are mangled, as in the Cartesian output, so reach them through an extern "C" shim)."""
    import ctypes
    import subprocess
    rep = generate_sh_code(5, "shcpp", language="c++", include_dir=str(tmp_path), src_dir=str(tmp_path))
    assert rep["source"].endswith(".cpp")
    shim = tmp_path / "shim.cpp"
    shim.write_text('#include "shcpp.h"\nextern "C" void shim_m2l(double x, double y, double z, double *M, double *L)'
                    " { M2L_4(x, y, z, M, L); }\n")
    so = str(tmp_path / "shcpp.so")
    subprocess.run(["g++", "-O1", "-Wall", "-shared", "-fPIC", "-o", so, rep["source"], str(shim), "-lm"], check=True)
    lib = ctypes.CDLL(so)
    dp = ctypes.POINTER(ctypes.c_double)
    lib.shim_m2l.argtypes = [ctypes.c_double] * 3 + [dp, dp]
    rng = np.random.default_rng(32)
    M = rng.standard_normal(25)
    L = np.zeros(25)
    lib.shim_m2l(2.0, 1.0, -2.5, M.ctypes.data_as(dp), L.ctypes.data_as(dp))
    ref = np.zeros(25)
    sm.m2l_rotation(4, 2.0, 1.0, -2.5, M, ref)
    assert rel_err(L, ref) < 1e-12


def test_header_sizes(tmp_path):
    rep = generate_sh_code(7, "hs", include_dir=str(tmp_path))
    h = open(rep["header"]).read()
    assert "FMMGEN_MULTIPOLESIZE[] = {1, 4, 9, 16, 25, 36, 49, 64}" in h
    assert "FMMGEN_MAXORDER 7" in h and "M2L_6(" in h and "M2L_7(" not in h and "S2M_1(" in h
