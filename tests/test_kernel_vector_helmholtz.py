"""Stokeslet API and low-frequency Helmholtz kernels."""
import shutil
import subprocess

import numpy as np
import pytest
import sympy as sp

from fmmgen import kernel as K
from fmmgen import Kernel
from fmmgen.kernel import nterms, source_indices
from fmmgen.vector import VectorKernel, kelvin, stokeslet
from kernel_util import HAVE_GCC, Lib, cluster, far_field, rel_err
from test_kernel_derivatives import evaluate_D, reference

PT = [0.7, -0.4, 0.9]


@pytest.mark.parametrize("make", [K.helmholtz_cos, K.helmholtz_sin])
def test_helmholtz_derivatives_and_pde(make):
    k = make(1.3)  # the constructor checks (Laplace + k^2) G = 0
    got, ref = evaluate_D(k, 4, PT), reference(k, 4, PT)
    for n in ref:
        assert got[n] == pytest.approx(ref[n], rel=1e-8, abs=1e-9), n


@pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")
@pytest.mark.parametrize("make", [K.helmholtz_cos, K.helmholtz_sin])
def test_helmholtz_far_field_converges(make, tmp_path):
    lib = Lib(make(1.0), tmp_path, [2, 4, 6, 8])
    errs = []
    for p in (2, 4, 6, 8):
        phi, F, pr, Fr = far_field(lib, p, np.random.default_rng(3))
        errs.append(rel_err(phi, pr))
    assert errs[3] < 1e-3 * errs[0], errs


def _stokes_array():
    x, y, z = X = sp.symbols("x y z")
    r = sp.sqrt(x**2 + y**2 + z**2)
    return np.eye(3, dtype=int) / r + np.einsum("i,j->ij", np.array(X, dtype=object), np.array(X, dtype=object)) / r**3


def test_kernel_dispatches_on_arrays():
    k = Kernel(_stokes_array())
    assert isinstance(k, VectorKernel) and (k.n_out, k.n_src) == (3, 3) and k.contract == "ij,j->i"
    assert isinstance(Kernel(1 / sp.Symbol("x")), Kernel)
    # split: 3 diagonal entries + 3 symmetric off-diagonal pairs
    assert len(k.groups) == 6 and sum(g.nchannels for g in k.groups) == 9
    bad = pytest.raises(ValueError)
    with bad:
        Kernel(_stokes_array(), contract="ij,k->i")


def test_decomposition_matches_array():
    ref = Kernel(_stokes_array())
    st = stokeslet()
    rng = np.random.default_rng(0)
    src, tgt, f = rng.standard_normal((6, 3)), rng.standard_normal((3, 3)) + 4, rng.standard_normal((6, 3))
    assert np.allclose(st.direct(tgt, src, f), ref.direct(tgt, src, f), rtol=1e-12)
    assert len(st.groups) == 1 and st.groups[0].nchannels == 3
    with pytest.raises(ValueError):  # a wrong decomposition is rejected at construction
        VectorKernel(_stokes_array(), decomposition=[kelvin(0.25).groups[0]], split_channels=False)


def _far_field(vk, p, tmp_path, rng):
    d = vk.dim
    cs, ct = np.zeros(3), np.array([2.5, 0.5, -0.5])
    src = cluster(rng, cs, 10, 0.3, 3)
    f = rng.standard_normal((10, vk.n_src))
    tgt = cluster(rng, ct, 4, 0.3, 3)
    u = np.zeros((4, vk.n_out))
    for g, grp in enumerate(vk.groups):
        (tmp_path / f"g{g}").mkdir(parents=True, exist_ok=True)
        lib = Lib(grp.kernel, tmp_path / f"g{g}", [p], source_order=grp.source_order, nchannels=grp.nchannels)
        S = vk.source_weights(g, f)
        M = np.zeros(len(source_indices(d, p, grp.source_order)) * grp.nchannels)
        for y_, Sj in zip(src, S):
            lib.call(f"S2M_{p}", cs - y_, np.ascontiguousarray(Sj), M)
        L = np.zeros(nterms(d, p - grp.source_order) * grp.nchannels)
        lib.call(f"M2L_{p}", ct - cs, M, L)
        for i, t in enumerate(tgt):
            F = np.zeros(grp.nchannels * (1 + d))
            lib.call(f"L2P_{p}", t - ct, L, F)
            vk.output(g, F, u[i])
    return u, vk.direct(tgt, src, f)


@pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")
@pytest.mark.parametrize("name", ["decomposition", "split"])
def test_stokeslet_far_field(name, tmp_path):
    vk = stokeslet() if name == "decomposition" else Kernel(_stokes_array(), method="symbolic")
    assert (len(vk.groups) == 1) == (name == "decomposition")
    e4 = rel_err(*_far_field(vk, 4, tmp_path / "p4", np.random.default_rng(2)))
    e6 = rel_err(*_far_field(vk, 6, tmp_path / "p6", np.random.default_rng(2)))
    assert e6 < 0.5 * e4 and e6 < 2e-2, (name, e4, e6)


def _split_stokeslet():
    k = Kernel(_stokes_array(), method="symbolic")
    assert len(k.groups) == 6
    return k


@pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")
def test_groups_link_into_one_program(tmp_path):
    # every group defines S2M_p, P2P, ...: without per-group prefixes they collide at link time
    _split_stokeslet().generate("st", 2, orders=[2], src_dir=str(tmp_path), include_dir=str(tmp_path))
    cfiles = sorted(str(p) for p in tmp_path.glob("st_g*.c"))
    assert len(cfiles) == 6 and "st_g3_M2L_2" in (tmp_path / "st_g3.h").read_text()
    subprocess.run(["gcc", "-O0", "-shared", "-fPIC", "-I", str(tmp_path), *cfiles, "-o", str(tmp_path / "st.so"),
                    "-lm"], check=True, capture_output=True)


@pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")
def test_fortran_groups_are_separate_modules(tmp_path):
    _split_stokeslet().generate("st", 2, orders=[2], language="fortran", src_dir=str(tmp_path))
    files = sorted(p.name for p in tmp_path.glob("st_g*.f90"))
    assert len(files) == 6
    subprocess.run(["gfortran", "-std=f2008", "-c", *files], cwd=tmp_path, check=True, capture_output=True)
