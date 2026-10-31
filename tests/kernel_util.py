"""Shared helpers for the test_kernel_*.py tests: generate, compile (gcc, as a
shared library) and drive generated operators from Python through ctypes."""
import ctypes
import os
import shutil
import subprocess

import numpy as np
import sympy as sp

from fmmgen.kernel import Kernel, generate_kernel_code, nterms

x, y, z = sp.symbols("x y z")
r3 = sp.sqrt(x**2 + y**2 + z**2)
r2 = sp.sqrt(x**2 + y**2)

LAP3 = {(2, 0, 0): 1, (0, 2, 0): 1, (0, 0, 2): 1}
LAP2 = {(2, 0): 1, (0, 2): 1}
BIH3 = {(4, 0, 0): 1, (0, 4, 0): 1, (0, 0, 4): 1, (2, 2, 0): 2, (2, 0, 2): 2, (0, 2, 2): 2}

# name -> factory. Each returns a fresh Kernel.
KERNELS = {
    "laplace": lambda **kw: Kernel(1 / r3, pde=LAP3, **kw),
    "yukawa": lambda **kw: Kernel(
        sp.exp(-sp.Symbol("k") * r3) / r3, params={"k": 0.7},
        pde={**LAP3, (0, 0, 0): -sp.Symbol("k") ** 2}, **kw),
    "log2d": lambda **kw: Kernel(sp.log(r2), (x, y), pde=LAP2, **kw),
    "r": lambda **kw: Kernel(r3, pde=BIH3, **kw),
    "dipolelike": lambda **kw: Kernel(1 / r3**3, **kw),
    "plummer": lambda **kw: Kernel(1 / sp.sqrt(x**2 + y**2 + z**2 + sp.Symbol("eps") ** 2),
                                   params={"eps": 0.3}, **kw),
    "gaussian": lambda **kw: Kernel(sp.exp(-(x**2 + y**2 + z**2)), **kw),
    "aniso": lambda **kw: Kernel(1 / sp.sqrt(x**2 + 2 * y**2 + 3 * z**2), **kw),
}

HAVE_GCC = shutil.which("gcc") is not None
HAVE_GXX = shutil.which("g++") is not None


class Lib:
    """Compiled operators for one kernel."""

    def __init__(self, kernel, tmpdir, orders, language="c", **gen):
        self.kernel = kernel
        self.dim = kernel.dim
        self.report = generate_kernel_code(kernel, max(orders), "ops", orders=orders, language=language,
                                           include_dir=str(tmpdir), src_dir=str(tmpdir), **gen)
        src = os.path.join(str(tmpdir), "ops." + ("c" if language == "c" else "cpp"))
        so = os.path.join(str(tmpdir), "ops.so")
        cc = "gcc" if language == "c" else "g++"
        subprocess.run([cc, "-O1", "-shared", "-fPIC", "-o", so, src, "-lm"], check=True)
        self.lib = ctypes.CDLL(so)

    def call(self, fname, coords, *arrays):
        f = getattr(self.lib, fname)
        dp = ctypes.POINTER(ctypes.c_double)
        f.argtypes = [ctypes.c_double] * self.dim + [dp] * len(arrays)
        f.restype = None
        f(*[float(c) for c in coords], *[a.ctypes.data_as(dp) for a in arrays])


def direct(kernel, targets, sources, q):
    """phi and field at targets by direct summation (numpy)."""
    G, grads = kernel.lambdify()
    phi = np.zeros(len(targets))
    F = np.zeros((len(targets), kernel.dim))
    for i, t in enumerate(targets):
        d = t[None, :] - sources
        cols = [d[:, j] for j in range(kernel.dim)]
        phi[i] = np.sum(q * G(*cols))
        for j in range(kernel.dim):
            F[i, j] = -np.sum(q * grads[j](*cols))
    return phi, F


def cluster(rng, centre, n, radius, dim):
    pts = rng.standard_normal((n, dim))
    pts *= (radius * rng.random((n, 1)) ** (1 / dim)) / np.linalg.norm(pts, axis=1, keepdims=True)
    return centre + pts


def far_field(lib, p, rng, route="A", nsrc=20, ntgt=8, dist=3.0, radius=0.4):
    """Run P2M -> [M2M] -> M2L -> [L2L] -> L2P; return (phi, F, phi_ref, F_ref)."""
    d = lib.dim
    cs = rng.standard_normal(d) * 0.1
    u = rng.standard_normal(d)
    u /= np.linalg.norm(u)
    ct = cs + dist * u
    N = nterms(d, p)
    osz = 1 + d
    q = rng.standard_normal(nsrc)
    if route == "A":
        src = cluster(rng, cs, nsrc, radius, d)
        tgt = cluster(rng, ct, ntgt, radius, d)
        M = np.zeros(N)
        for s, qq in zip(src, q):
            lib.call(f"S2M_{p}", cs - s, np.array([qq]), M)
        L = np.zeros(N)
        lib.call(f"M2L_{p}", ct - cs, M, L)
        F = np.zeros((ntgt, osz))
        for i, t in enumerate(tgt):
            lib.call(f"L2P_{p}", t - ct, L, F[i])
    else:
        # two child source cells -> M2M -> M2L -> L2L to two child target cells
        off = radius / 2 * np.ones(d) / np.sqrt(d)
        cs_kids = [cs + off, cs - off]
        ct_kids = [ct + off, ct - off]
        half = radius / 2
        src = np.vstack([cluster(rng, c, nsrc // 2, half, d) for c in cs_kids])
        q = q[: len(src)]
        M = np.zeros(N)
        for k, c in enumerate(cs_kids):
            Mk = np.zeros(N)
            for s, qq in zip(src[k * (nsrc // 2):(k + 1) * (nsrc // 2)], q[k * (nsrc // 2):(k + 1) * (nsrc // 2)]):
                lib.call(f"S2M_{p}", c - s, np.array([qq]), Mk)
            lib.call(f"M2M_{p}", cs - c, Mk, M)
        L = np.zeros(N)
        lib.call(f"M2L_{p}", ct - cs, M, L)
        tgt = np.vstack([cluster(rng, c, ntgt // 2, half, d) for c in ct_kids])
        F = np.zeros((len(tgt), osz))
        for k, c in enumerate(ct_kids):
            Lk = np.zeros(N)
            lib.call(f"L2L_{p}", c - ct, L, Lk)
            for i in range(k * (ntgt // 2), (k + 1) * (ntgt // 2)):
                lib.call(f"L2P_{p}", tgt[i] - c, Lk, F[i])
    phi_ref, F_ref = direct(lib.kernel, tgt, src, q)
    return F[:, 0], F[:, 1:], phi_ref, F_ref


def rel_err(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


# ---------------------------------------------------------------------------
# Multipole sources (source_order s >= 0), possibly several channels
# ---------------------------------------------------------------------------
from fmmgen.kernel import source_indices, shell_indices  # noqa: E402


def source_derivs(kernel, s, field=True):
    """Lambdified d^m G for |m| = s (and m + e_i for the field)."""
    d = kernel.dim
    sh = shell_indices(d, s)
    pot, fld = {}, {}

    def dd(m):
        e = kernel.expr
        for v, k in zip(kernel.variables, m):
            if k:
                e = sp.diff(e, v, k)
        return kernel.lambdify_expr(e)

    for m in sh:
        pot[m] = dd(m)
        if field:
            for i in range(d):
                mi = tuple(k + (j == i) for j, k in enumerate(m))
                fld[(m, i)] = dd(mi)
    return sh, pot, fld


def direct_s(kernel, targets, sources, S, s):
    """phi, F at targets for sources with coefficients S[j, shell] (order s)."""
    sh, pot, fld = source_derivs(kernel, s)
    d = kernel.dim
    phi = np.zeros(len(targets))
    F = np.zeros((len(targets), d))
    for i, t in enumerate(targets):
        dd = t[None, :] - sources
        cols = [dd[:, j] for j in range(d)]
        for m, im in sh.items():
            phi[i] += np.sum(S[:, im] * pot[m](*cols))
            for a in range(d):
                F[i, a] -= np.sum(S[:, im] * fld[(m, a)](*cols))
    return phi, F


def far_field_s(lib, p, s, rng, dist=3.0, radius=0.4, nsrc=12, ntgt=6):
    """S2M -> M2L -> L2P for order-s sources; returns (phi, F, phi_ref, F_ref)."""
    d = lib.dim
    sh = shell_indices(d, s)
    NM = len(source_indices(d, p, s))
    NL = nterms(d, p - s)
    cs = rng.standard_normal(d) * 0.1
    u = rng.standard_normal(d)
    u /= np.linalg.norm(u)
    ct = cs + dist * u
    src = cluster(rng, cs, nsrc, radius, d)
    S = rng.standard_normal((nsrc, len(sh)))
    tgt = cluster(rng, ct, ntgt, radius, d)
    M = np.zeros(NM)
    for y_, Sj in zip(src, S):
        lib.call(f"S2M_{p}", cs - y_, np.ascontiguousarray(Sj), M)
    L = np.zeros(NL)
    lib.call(f"M2L_{p}", ct - cs, M, L)
    F = np.zeros((ntgt, 1 + d))
    for i, t in enumerate(tgt):
        lib.call(f"L2P_{p}", t - ct, L, F[i])
    phi_ref, F_ref = direct_s(lib.kernel, tgt, src, S, s)
    return F[:, 0], F[:, 1:], phi_ref, F_ref
