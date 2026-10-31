"""Separable multi-channel kernels: Lennard-Jones with geometric and
Lorentz-Berthelot mixing, per-particle (eps, sigma), compiled operators
against direct summation."""
import shutil
import subprocess

import numpy as np
import pytest

from fmmgen.kernel import nterms, source_indices
from fmmgen.separable import SeparableKernel
from kernel_util import HAVE_GCC, Lib, cluster

pytestmark = pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")


def far_field_separable(sk, libs, p, rng, dist=2.6, radius=0.4, nsrc=14, ntgt=5):
    d = 3
    cs = np.zeros(d)
    ct = np.array([1.0, 0.4, -0.3]) / np.linalg.norm([1.0, 0.4, -0.3]) * dist
    src, tgt = cluster(rng, cs, nsrc, radius, d), cluster(rng, ct, ntgt, radius, d)
    th_s = np.column_stack([rng.uniform(0.5, 1.5, nsrc), rng.uniform(0.8, 1.2, nsrc)])
    th_t = np.column_stack([rng.uniform(0.5, 1.5, ntgt), rng.uniform(0.8, 1.2, ntgt)])
    phi = np.zeros(ntgt)
    F = np.zeros((ntgt, d))
    for g, lib in enumerate(libs):
        C = len(sk.groups[g][1])
        NM, NL = len(source_indices(d, p, 0)), nterms(d, p)
        b = sk.source_weights(g, th_s)
        M = np.zeros(NM * C)
        for y_, bj in zip(src, b):
            lib.call(f"S2M_{p}", cs - y_, np.ascontiguousarray(bj), M)
        L = np.zeros(NL * C)
        lib.call(f"M2L_{p}", ct - cs, M, L)
        a = sk.target_weights(g, th_t)
        for i, t in enumerate(tgt):
            out = np.zeros(4 * C)
            lib.call(f"L2P_{p}", t - ct, L, out)
            out = out.reshape(C, 4)
            phi[i] += np.sum(a[i] * out[:, 0])
            F[i] += (a[i][:, None] * out[:, 1:]).sum(axis=0)
    phi_ref, F_ref = sk.direct(tgt, th_t, src, th_s)
    return phi, F, phi_ref, F_ref


@pytest.mark.parametrize("mixing,nch", [("geometric", 2), ("lorentz-berthelot", 20)])
def test_lennard_jones_mixing(mixing, nch, tmp_path):
    sk = SeparableKernel.lennard_jones(mixing)
    assert sk.nchannels == nch
    assert len(sk.groups) == 2  # one derivative tensor per exponent, shared across channels
    libs = []
    for g, (k, ts) in enumerate(sk.groups):
        sub = tmp_path / f"g{g}"
        sub.mkdir()
        libs.append(Lib(k, sub, [3, 6], nchannels=len(ts)))
    errs = {}
    for p in (3, 6):
        phi, F, pr, Fr = far_field_separable(sk, libs, p, np.random.default_rng(2))
        errs[p] = (np.max(np.abs(phi - pr)) / np.max(np.abs(pr)), np.max(np.abs(F - Fr)) / np.max(np.abs(Fr)))
    assert errs[6][0] < 0.1 * errs[3][0] and errs[6][1] < 0.1 * errs[3][1], errs
    assert errs[6][0] < 5e-3, errs


def _stokes_S(f, shell):
    """Per-output-component source coefficients for the Stokeslet in terms of
    B = r: u_i = [f_i Lap - sum_j f_j d_i d_j] B / (8 pi mu) (mu = 1/(8 pi))."""
    S = np.zeros((3, len(shell)))
    for i in range(3):
        for m, im in shell.items():
            if max(m) == 2:  # d_k^2
                k = m.index(2)
                S[i, im] = f[i] if k != i else 0.0
            else:  # d_a d_b, a != b
                a, b = [q for q, e in enumerate(m) if e]
                if a == i:
                    S[i, im] = -f[b]
                elif b == i:
                    S[i, im] = -f[a]
    return S


def _stokeslet(x, f):
    r = np.linalg.norm(x)
    return f / r + x * np.dot(x, f) / r**3  # mu = 1/(8 pi)


def test_stokeslet_from_biharmonic_kernel(tmp_path):
    """A vector-valued kernel from the existing machinery: the Stokeslet is
    three order-2 channels of the scalar kernel r (shared derivative tensor)."""
    from fmmgen.kernel import Kernel, shell_indices
    from kernel_util import BIH3, r3

    k = Kernel(r3, pde=BIH3)
    shell = shell_indices(3, 2)
    p = 6
    lib = Lib(k, tmp_path, [p], source_order=2, nchannels=3)
    rng = np.random.default_rng(1)
    src = cluster(rng, np.zeros(3), 8, 0.3, 3)
    f = rng.standard_normal((8, 3))
    ct = np.array([2.5, 0.5, -0.5])
    tgt = cluster(rng, ct, 3, 0.3, 3)
    t = tgt[0]
    u = np.zeros(3)
    for y_, fj in zip(src, f):
        F = np.zeros(12)
        lib.call("P2P", t - y_, np.ascontiguousarray(_stokes_S(fj, shell).reshape(-1)), F)
        u += F.reshape(3, 4)[:, 0]
    ref = sum(_stokeslet(t - y_, fj) for y_, fj in zip(src, f))
    assert np.allclose(u, ref, rtol=1e-9)
    NM = len(source_indices(3, p, 2))
    NL = nterms(3, p - 2)
    M = np.zeros(NM * 3)
    for y_, fj in zip(src, f):
        lib.call(f"S2M_{p}", -y_, np.ascontiguousarray(_stokes_S(fj, shell).reshape(-1)), M)
    L = np.zeros(NL * 3)
    lib.call(f"M2L_{p}", ct, M, L)
    for t in tgt:
        F = np.zeros(12)
        lib.call(f"L2P_{p}", t - ct, L, F)
        ref = sum(_stokeslet(t - y_, fj) for y_, fj in zip(src, f))
        assert np.allclose(F.reshape(3, 4)[:, 0], ref, rtol=2e-2, atol=1e-3 * np.max(np.abs(ref)))


@pytest.mark.parametrize("method", ["svd", "chebyshev"])
@pytest.mark.parametrize("n", [6, 12])
def test_low_rank_sigma_factors(method, n):
    from fmmgen.separable import low_rank_sigma_factors

    rng = np.random.default_rng(0)
    sig = rng.uniform(2.5, 4.0, 300)
    exact = ((sig[:, None] + sig[None, :]) / 2) ** n
    prev = None
    for tol in (1e-2, 1e-4, 1e-6):
        fa, gb, rank, err = low_rank_sigma_factors(n, sigma_range=(2.5, 4.0), tol=tol, method=method)
        approx = fa(sig) @ gb(sig).T
        assert np.max(np.abs(approx - exact) / exact) < 2 * tol  # holds off-grid too
        assert prev is None or rank >= prev
        prev = rank
    assert rank <= n + 1
    # rank n+1 reproduces the exact binomial expansion to rounding
    fa, gb, _, err = low_rank_sigma_factors(n, sigma_range=(2.5, 4.0), rank=n + 1, method="chebyshev")
    assert err < 1e-9


def test_low_rank_lj_weights_match_exact_pair_coefficient():
    """Sum over channels of a_c(i) b_c(j) reproduces the exact LB pair coefficient."""
    exact = SeparableKernel.lennard_jones("lorentz-berthelot")
    comp = SeparableKernel.lennard_jones_lowrank(sigma_range=(2.5, 4.0), tol=1e-5)
    assert comp.nchannels < exact.nchannels
    rng = np.random.default_rng(1)
    th = np.column_stack([rng.uniform(0.02, 0.25, 40), rng.uniform(2.5, 4.0, 40)])
    for g in (0, 1):
        # same exponent groups in the same order (12 then 6)
        a_e, b_e = exact.target_weights(g, th), exact.source_weights(g, th)
        a_c, b_c = comp.target_weights(g, th), comp.source_weights(g, th)
        full_e = a_e @ b_e.T
        full_c = a_c @ b_c.T
        assert np.max(np.abs(full_c - full_e) / np.abs(full_e)) < 1e-4


def test_groups_link_into_one_program(tmp_path):
    # the r^-12 and r^-6 groups both define S2M_p, P2P, ...: they need their own namespaces
    sk = SeparableKernel.lennard_jones("geometric")
    sk.generate("lj", 2, orders=[2], src_dir=str(tmp_path), include_dir=str(tmp_path))
    cfiles = sorted(str(p) for p in tmp_path.glob("lj_g*.c"))
    assert len(cfiles) == len(sk.groups) == 2
    subprocess.run(["gcc", "-O0", "-shared", "-fPIC", "-I", str(tmp_path), *cfiles, "-o", str(tmp_path / "lj.so"),
                    "-lm"], check=True, capture_output=True)


@pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not available")
def test_fortran_groups_are_separate_modules(tmp_path):
    sk = SeparableKernel.lennard_jones("geometric")
    sk.generate("lj", 2, orders=[2], language="fortran", src_dir=str(tmp_path))
    files = sorted(p.name for p in tmp_path.glob("lj_g*.f90"))
    assert len(files) == 2
    subprocess.run(["gfortran", "-std=f2008", "-c", *files], cwd=tmp_path, check=True, capture_output=True)


@pytest.mark.parametrize("nsamples", [60, 900])  # 900 distinct values takes the quantile-grid branch
def test_low_rank_from_sampled_sigma_and_fixed_ranks(nsamples):
    from fmmgen.separable import low_rank_sigma_factors

    sig = np.random.default_rng(3).uniform(2.5, 4.0, nsamples)
    exact = ((sig[:, None] + sig[None, :]) / 2) ** 6
    fa, gb, rank, err = low_rank_sigma_factors(6, sigma_sample=sig, tol=1e-5, method="svd")
    assert err < 1e-5 and rank < 7
    assert np.max(np.abs(fa(sig) @ gb(sig).T - exact) / exact) < 1e-4  # the samples themselves
    sk = SeparableKernel.lennard_jones_lowrank(sigma_sample=sig, ranks=(5, 3))
    assert sk.compression[12][0] == 5 and sk.compression[6][0] == 3 and sk.nchannels == 8


def test_bad_mixing_and_method_are_rejected():
    from fmmgen.separable import low_rank_sigma_factors

    with pytest.raises(ValueError, match="mixing"):
        SeparableKernel.lennard_jones("arithmetic")
    with pytest.raises(ValueError, match="method"):
        low_rank_sigma_factors(6, sigma_range=(2.5, 4.0), method="qr")
