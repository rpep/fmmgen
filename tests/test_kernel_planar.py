"""Planar variants (sources and targets in z = 0) of the kernel-generic operators: on in-plane points they must
equal the full operators restricted to the retained array entries, for any kernel."""
import shutil
import subprocess

import numpy as np
import pytest
import sympy as sp

from fmmgen.kernel import (Kernel, _planar_keep, generate_kernel_code, nterms, planar_sizes, power_law,
                           shell_indices, source_indices, yukawa)
from kernel_util import HAVE_GCC, Lib

x, y, z = sp.symbols("x y z")
r = sp.sqrt(x**2 + y**2 + z**2)

CASES = [
    ("yukawa-charges", lambda: yukawa(1.0), 0, 1, 3),
    ("r^-3-dipoles", lambda: power_law(3), 1, 1, 3),
    ("yukawa-quadrupoles-2ch", lambda: yukawa(0.7), 2, 2, 4),
    # not even in z, not radial: nothing may assume the symmetry of the Laplace kernel
    ("odd-in-z", lambda: Kernel(sp.exp(-r) * (1 + z / (r + 1)) / r, (x, y, z)), 0, 1, 3),
]


@pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")
@pytest.mark.parametrize("name,make,s,nch,p", CASES, ids=[c[0] for c in CASES])
def test_planar_equals_full_restricted_to_the_plane(tmp_path, name, make, s, nch, p):
    d = 3
    lib = Lib(make(), tmp_path, [p], source_order=s, nchannels=nch, planar=True)
    km, kl = _planar_keep(d, p, s, nch, True)
    nm, nl = len(source_indices(d, p, s)) * nch, nterms(d, p - s) * nch
    assert (len(km), len(kl)) == planar_sizes(d, p, s, nch, True) and len(km) < nm and len(kl) < nl
    rng = np.random.default_rng(7)
    shell, nf = len(shell_indices(d, s)) * nch, 4 * nch
    S = rng.standard_normal(shell)

    def same(full, xy, rows=None, what=""):
        ref = full if rows is None else full[rows]
        assert np.allclose(xy, ref, rtol=1e-11, atol=1e-11 * np.max(np.abs(full))), what

    c = [(0.3, -0.2, 0.0), (-0.5, 0.4, 0.0), (2.0, -1.5, 0.0), (0.4, 0.1, 0.0), (0.2, -0.1, 0.0)]
    # S2M, and the entries the planar array drops really are zero
    M = np.zeros(nm)
    lib.call(f"S2M_{p}", c[0], S, M)
    Mxy = np.zeros(len(km))
    lib.call(f"S2Mxy_{p}", c[0], S, Mxy)
    same(M, Mxy, km, "S2M")
    assert not np.any(np.delete(M, km)), "S2M of in-plane sources has entries outside the planar set"
    # M2M
    Ms, Msxy = np.zeros(nm), np.zeros(len(km))
    lib.call(f"M2M_{p}", c[1], M, Ms)
    lib.call(f"M2Mxy_{p}", c[1], Mxy, Msxy)
    same(Ms, Msxy, km, "M2M")
    # M2L
    L, Lxy = np.zeros(nl), np.zeros(len(kl))
    lib.call(f"M2L_{p}", c[2], Ms, L)
    lib.call(f"M2Lxy_{p}", c[2], Msxy, Lxy)
    same(L, Lxy, kl, "M2L")
    # L2L, L2P (the full local array has non-zero entries outside the planar set: they must not matter)
    Ls, Lsxy = np.zeros(nl), np.zeros(len(kl))
    lib.call(f"L2L_{p}", c[3], L, Ls)
    lib.call(f"L2Lxy_{p}", c[3], L[kl], Lsxy)
    same(Ls, Lsxy, kl, "L2L")
    F, Fxy = np.zeros(nf), np.zeros(nf)
    lib.call(f"L2P_{p}", c[4], L, F)
    lib.call(f"L2Pxy_{p}", c[4], L[kl], Fxy)
    same(F, Fxy, None, "L2P")
    # M2P, P2P
    F, Fxy = np.zeros(nf), np.zeros(nf)
    lib.call(f"M2P_{p}", c[2], M, F)
    lib.call(f"M2Pxy_{p}", c[2], Mxy, Fxy)
    same(F, Fxy, None, "M2P")
    F, Fxy = np.zeros(nf), np.zeros(nf)
    lib.call("P2P", c[2], S, F)
    lib.call("P2Pxy", c[2], S, Fxy)
    same(F, Fxy, None, "P2P")


def test_planar_header_tables_and_wrappers(tmp_path):
    generate_kernel_code(yukawa(1.0), 4, "k", orders=[2, 4], planar=True, prefix="y_", src_dir=str(tmp_path),
                         include_dir=str(tmp_path))
    h = (tmp_path / "k.h").read_text()
    assert "#define FMMGEN_Y_PLANAR 1" in h
    # the tables are indexed by order 0..order
    want_m = ", ".join(str(planar_sizes(3, o, 0, 1, True)[0]) for o in range(5))
    want_l = ", ".join(str(planar_sizes(3, o, 0, 1, True)[1]) for o in range(5))
    assert f"FMMGEN_Y_PLANAR_MULTIPOLESIZE[] = {{{want_m}}}" in h
    assert f"FMMGEN_Y_PLANAR_LOCALSIZE[] = {{{want_l}}}" in h
    assert "void y_M2Lxy_4(" in h and "void y_M2Lxy(" in h and "void y_M2Lxy_3(" not in h
    assert "void y_P2Pxy(" in h and "void y_P2P_batchxy(" in h
    with pytest.raises(ValueError, match="three-dimensional"):
        generate_kernel_code(power_law(1, dim=2), 2, "k2", planar=True, src_dir=str(tmp_path),
                             include_dir=str(tmp_path))


@pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")
def test_planar_fortran_compiles(tmp_path):
    generate_kernel_code(yukawa(1.0), 3, "kp", orders=[3], planar=True, language="fortran", src_dir=str(tmp_path))
    text = (tmp_path / "kp.f90").read_text()
    assert "FMMGEN_PLANAR = .true." in text and "subroutine M2Lxy_3" in text and "subroutine M2Lxy(" in text
    subprocess.run(["gfortran", "-std=f2008", "-fsyntax-only", "kp.f90"], cwd=tmp_path, check=True,
                   capture_output=True)
