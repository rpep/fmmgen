"""Helpers for the spherical-harmonic operator tests: generate, compile (gcc, shared
library) and drive the generated operators through ctypes, in the same style as
kernel_util.Lib."""
import ctypes
import os
import subprocess

import numpy as np

from fmmgen.sh import generate_sh_code
from fmmgen import sh_math as sm


class SHLib:
    def __init__(self, tmpdir, order, name="sh", opt="-O1", language="c", **gen):
        self.order = order
        self.report = generate_sh_code(order, name, include_dir=str(tmpdir), src_dir=str(tmpdir),
                                       language=language, **gen)
        so = os.path.join(str(tmpdir), name + ".so")
        cc = "gcc" if language == "c" else "g++"
        subprocess.run([cc, opt, "-shared", "-fPIC", "-o", so, self.report["source"], "-lm"], check=True)
        self.lib = ctypes.CDLL(so)

    def call(self, fname, coords, *arrays):
        f = getattr(self.lib, fname)
        dp = ctypes.POINTER(ctypes.c_double)
        f.argtypes = [ctypes.c_double] * len(coords) + [dp] * len(arrays)
        f.restype = None
        f(*[float(c) for c in coords], *[np.ascontiguousarray(a).ctypes.data_as(dp) for a in arrays])

    def apply(self, fname, *arrays):
        f = getattr(self.lib, fname)
        dp = ctypes.POINTER(ctypes.c_double)
        f.argtypes = [dp] * len(arrays)
        f.restype = None
        f(*[a.ctypes.data_as(dp) for a in arrays])


def direct(targets, sources, q):
    """phi and F = -grad phi at targets for 1/r charges (numpy)."""
    phi = np.zeros(len(targets))
    F = np.zeros((len(targets), 3))
    for i, t in enumerate(targets):
        d = t[None, :] - sources
        r = np.linalg.norm(d, axis=1)
        phi[i] = np.sum(q / r)
        F[i] = np.sum(q[:, None] * d / r[:, None] ** 3, axis=0)
    return phi, F


def cluster(rng, centre, n, half):
    return centre + rng.uniform(-half, half, (n, 3))


def far_field_chain(lib, p, rng, m2l="M2L", nsrc=20, ntgt=6, dist=(3.0, 1.0, -2.0), half=0.5, route="A"):
    """S2M -> [M2M] -> M2L -> [L2L] -> L2P using the compiled library; returns
    (phi, F, phi_ref, F_ref)."""
    N = sm.sizes(p)[0]
    cs = np.zeros(3)
    ct = np.array(dist, dtype=float)
    q = rng.standard_normal(nsrc)
    if route == "A":
        src = cluster(rng, cs, nsrc, half)
        tgt = cluster(rng, ct, ntgt, half)
        M = np.zeros(N)
        for s, qq in zip(src, q):
            lib.call(f"S2M_{p}", cs - s, np.array([qq]), M)
        L = np.zeros(N)
        lib.call(f"{m2l}_{p}", ct - cs, M, L)
        F = np.zeros((ntgt, 4))
        for i, t in enumerate(tgt):
            lib.call(f"L2P_{p}", t - ct, L, F[i])
    else:
        off = np.full(3, half / 2)
        cs_k = [cs + off, cs - off]
        ct_k = [ct + off, ct - off]
        src = np.vstack([cluster(rng, c, nsrc // 2, half / 2) for c in cs_k])
        q = q[: len(src)]
        M = np.zeros(N)
        for k, c in enumerate(cs_k):
            Mk = np.zeros(N)
            for s, qq in zip(src[k * (nsrc // 2):(k + 1) * (nsrc // 2)], q[k * (nsrc // 2):(k + 1) * (nsrc // 2)]):
                lib.call(f"S2M_{p}", c - s, np.array([qq]), Mk)
            lib.call(f"M2M_{p}", cs - c, Mk, M)
        L = np.zeros(N)
        lib.call(f"{m2l}_{p}", ct - cs, M, L)
        tgt = np.vstack([cluster(rng, c, ntgt // 2, half / 2) for c in ct_k])
        F = np.zeros((len(tgt), 4))
        for k, c in enumerate(ct_k):
            Lk = np.zeros(N)
            lib.call(f"L2L_{p}", c - ct, L, Lk)
            for i in range(k * (ntgt // 2), (k + 1) * (ntgt // 2)):
                lib.call(f"L2P_{p}", tgt[i] - c, Lk, F[i])
    phi_ref, F_ref = direct(tgt, src, q)
    return F[:, 0], F[:, 1:], phi_ref, F_ref


def rel_err(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))
