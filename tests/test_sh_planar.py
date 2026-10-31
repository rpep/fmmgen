"""Planar spherical-harmonic operators (charges, targets and centres in z = 0): they must equal the full
operators restricted to the n + |m| even entries, and converge to the direct sum."""
import re

import numpy as np
import pytest

from fmmgen import sh_math as sm
from fmmgen.sh import generate_sh_code
from sh_util import SHLib, cluster, direct, rel_err

pytestmark = pytest.mark.skipif(not __import__("kernel_util").HAVE_GCC, reason="gcc not available")


def test_planar_keep_is_the_even_parity_half():
    for p in range(6):
        keep = sm.planar_keep(p)
        assert len(keep) == sm.planar_size(p) == (p + 1) * (p + 2) // 2
        assert keep == sorted(set(keep)) and max(keep) < (p + 1) ** 2
        # n + |m| even, n + 1 per degree, and the planar index is n(n+1)/2 + (n+k)/2
        for pos, full in enumerate(keep):
            n = int(np.floor(np.sqrt(full)))
            k = full - n * n - n
            assert (n + abs(k)) % 2 == 0 and pos == n * (n + 1) // 2 + (n + k) // 2


@pytest.mark.parametrize("trunc,unroll", [("box", False), ("total", False), ("box", True)])
def test_planar_equals_full_restricted_to_the_plane(tmp_path, trunc, unroll):
    p = 5
    lib = SHLib(tmp_path, p + 1, planar=True, trunc=trunc, unroll=unroll)
    keep = sm.planar_keep(p)
    n, npl = sm.sizes(p)[0], sm.planar_size(p)
    c = [(0.3, -0.2, 0.0), (-0.5, 0.4, 0.0), (2.0, -1.5, 0.0), (0.4, 0.1, 0.0), (0.2, -0.1, 0.0)]

    def close(full, planar, what):
        assert np.allclose(planar, full, rtol=1e-11, atol=1e-11 * np.max(np.abs(full))), what

    q = np.array([0.7])
    M, Mxy = np.zeros(n), np.zeros(npl)
    lib.call(f"S2M_{p}", c[0], q, M)
    lib.call(f"S2Mxy_{p}", c[0], q, Mxy)
    close(M[keep], Mxy, "S2M")
    assert np.max(np.abs(np.delete(M, keep))) < 1e-13  # the entries a planar array drops are zero
    Ms, Msxy = np.zeros(n), np.zeros(npl)
    lib.call(f"M2M_{p}", c[1], M, Ms)
    lib.call(f"M2Mxy_{p}", c[1], Mxy, Msxy)
    close(Ms[keep], Msxy, "M2M")
    L, Lxy = np.zeros(n), np.zeros(npl)
    lib.call(f"M2L_{p}", c[2], Ms, L)
    lib.call(f"M2Lxy_{p}", c[2], Msxy, Lxy)
    close(L[keep], Lxy, "M2L")
    assert np.max(np.abs(np.delete(L, keep))) < 1e-12 * np.max(np.abs(L))
    Ls, Lsxy = np.zeros(n), np.zeros(npl)
    lib.call(f"L2L_{p}", c[3], L, Ls)
    lib.call(f"L2Lxy_{p}", c[3], Lxy, Lsxy)
    close(Ls[keep], Lsxy, "L2L")
    for fn_full, fn_xy, arr_full, arr_xy, xyz in ((f"L2P_{p}", f"L2Pxy_{p}", L, Lxy, c[4]),
                                                  (f"M2P_{p}", f"M2Pxy_{p}", Ms, Msxy, c[2])):
        F, Fxy = np.zeros(4), np.zeros(4)
        lib.call(fn_full, xyz, arr_full, F)
        lib.call(fn_xy, xyz, arr_xy, Fxy)
        close(F, Fxy, fn_xy)
        assert abs(Fxy[3]) < 1e-14 * np.max(np.abs(Fxy))  # no z field in the plane
    F, Fxy = np.zeros(4), np.zeros(4)
    lib.call("P2P", c[2], q, F)
    lib.call("P2Pxy", c[2], q, Fxy)
    close(F, Fxy, "P2P")


def _in_plane(rng, centre, n, half):
    pts = cluster(rng, centre, n, half)
    pts[:, 2] = 0.0
    return pts


def test_planar_fmm_chain_converges_to_the_direct_sum(tmp_path):
    lib = SHLib(tmp_path, 9, orders=[3, 8], planar=True)
    rng = np.random.default_rng(5)
    cs, ct, off = np.zeros(3), np.array([3.0, 1.0, 0.0]), np.array([0.25, 0.25, 0.0])
    src_cells, tgt_cells = (cs + off, cs - off), (ct + off, ct - off)  # two children of each parent cell
    errs = {}
    for p in (3, 8):
        npl = sm.planar_size(p)
        src = [_in_plane(rng, c, 10, 0.25) for c in src_cells]
        tgt = [_in_plane(rng, c, 3, 0.25) for c in tgt_cells]
        q = rng.standard_normal(20)
        M = np.zeros(npl)
        for k, c in enumerate(src_cells):
            Mk = np.zeros(npl)
            for s, qq in zip(src[k], q[10 * k:10 * k + 10]):
                lib.call(f"S2Mxy_{p}", c - s, np.array([qq]), Mk)
            lib.call(f"M2Mxy_{p}", cs - c, Mk, M)
        L = np.zeros(npl)
        lib.call(f"M2Lxy_{p}", ct - cs, M, L)
        F = np.zeros((6, 4))
        for k, c in enumerate(tgt_cells):
            Lk = np.zeros(npl)
            lib.call(f"L2Lxy_{p}", c - ct, L, Lk)
            for i, t in enumerate(tgt[k]):
                lib.call(f"L2Pxy_{p}", t - c, Lk, F[3 * k + i])
        phi_ref, F_ref = direct(np.vstack(tgt), np.vstack(src), q)
        errs[p] = (rel_err(F[:, 0], phi_ref), rel_err(F[:, 1:3], F_ref[:, :2]))
        assert np.all(F[:, 3] == 0.0)  # no z field in the plane
        assert np.allclose(F_ref[:, 2], 0.0)
    assert errs[8][0] < 1e-3 * errs[3][0] and errs[8][0] < 1e-6 and errs[8][1] < 1e-5, errs


def test_planar_header_tables_wrappers_and_prefix(tmp_path):
    generate_sh_code(7, "q", orders=[2, 5], planar=True, prefix="sp_", src_dir=str(tmp_path), include_dir=str(tmp_path))
    h = (tmp_path / "q.h").read_text()
    sizes = ", ".join(str(sm.planar_size(p)) for p in range(8))
    assert f"FMMGEN_SP_PLANAR_MULTIPOLESIZE[] = {{{sizes}}}" in h and "#define FMMGEN_SP_PLANAR 1" in h
    assert "void sp_M2Lxy_5(" in h and "void sp_M2Lxy(" in h and "void sp_P2Pxy(" in h and "sp_M2Lxy_3" not in h
    wrapper = re.search(r"void sp_M2Lxy\(.*?\n\}\n", (tmp_path / "q.c").read_text(), re.S).group(0)
    assert "case 2:" in wrapper and "case 5:" in wrapper and "case 3:" not in wrapper
    # without planar=True nothing planar is emitted
    generate_sh_code(4, "r", src_dir=str(tmp_path), include_dir=str(tmp_path))
    assert "xy" not in (tmp_path / "r.h").read_text() and "PLANAR" not in (tmp_path / "r.h").read_text()
