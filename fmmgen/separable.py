"""
Separable multi-channel kernels: U(i, j, r) = sum_c a_c(theta_i) b_c(theta_j) G_c(r).

Per-particle parameters theta (Lennard-Jones epsilon_i, sigma_i, charges, C6
coefficients, ...) enter a pair interaction through a mixing rule. Whenever the
rule makes the pair coefficient a sum of products of a function of i and a
function of j, the interaction is a short sum of ordinary radial kernels with
per-particle weights, and every channel is a plain FMM with "charges"
b_c(theta_j) and a target-side scale a_c(theta_i) applied after L2P:

    phi_i = sum_c a_c(theta_i) * [ sum_j b_c(theta_j) G_c(|x_i - x_j|) ].

* geometric mixing (sigma_ij = sqrt(sigma_i sigma_j), eps_ij = sqrt(eps_i eps_j)):
  Lennard-Jones is 2 channels, one r^-12 and one r^-6.
* Lorentz-Berthelot (sigma_ij = (sigma_i + sigma_j)/2): sigma_ij^n expands
  binomially into n+1 products, so LJ is 13 channels on r^-12 and 7 on r^-6.
  (n+1 is also the minimum: the tensor rank of (x+y)^n is n+1.)

Channels with the same G form a *group* and are generated as one multi-channel
operator set (`generate_kernel_code(..., nchannels=C)`): the derivative tensor
of the group is computed once and shared by the C contractions in M2L. That
saves the derivative work but not the contraction, which dominates.
"""
import math

import numpy as np
import sympy as sp

from .kernel import generate_kernel_code, power_law

__all__ = ["SeparableKernel", "low_rank_sigma_factors"]


def _cheb_nodes(lo, hi, m):
    k = np.arange(m)
    return 0.5 * (lo + hi) + 0.5 * (hi - lo) * np.cos((2 * k + 1) * np.pi / (2 * m))


def _lagrange(nodes, x):
    L = np.ones((len(x), len(nodes)))
    for k in range(len(nodes)):
        for j in range(len(nodes)):
            if j != k:
                L[:, k] *= (x - nodes[j]) / (nodes[k] - nodes[j])
    return L


def low_rank_sigma_factors(n, sigma_sample=None, sigma_range=None, tol=1e-4, rank=None, method="svd"):
    """Factors of ((x + y)/2)^n ~ sum_k fa(x)[:, k] gb(y)[:, k] over a bounded sigma range.

    Returns (fa, gb, rank, max_entrywise_relative_error_on_the_grid). fa and gb map an array of
    sigma to (len, rank) weight matrices.
    """
    if sigma_sample is not None:
        smp = np.asarray(sigma_sample, float)
        lo, hi = smp.min(), smp.max()
        grid = np.unique(smp) if len(np.unique(smp)) <= 400 else np.unique(np.quantile(smp, np.linspace(0, 1, 400)))
    else:
        lo, hi = sigma_range
        grid = np.linspace(lo, hi, 200)

    def pair(x, y):
        return ((np.asarray(x)[:, None] + np.asarray(y)[None, :]) / 2.0) ** n

    A = pair(grid, grid)
    if method == "chebyshev":
        def build(m):
            nodes = _cheb_nodes(lo, hi, m)
            return (lambda x: pair(x, nodes)), (lambda y: _lagrange(nodes, np.asarray(y, float)))
        maxr = n + 1
    elif method == "svd":
        lam, W = np.linalg.eigh(A)  # A is symmetric: A(x,y) = sum_k lam_k w_k(x) w_k(y)
        order = np.argsort(-np.abs(lam))
        lam, W = lam[order], W[:, order]

        def build(m):
            Wm, lm = W[:, :m], lam[:m]
            return (lambda x: pair(x, grid) @ Wm / lm), (lambda y: pair(y, grid) @ Wm)
        maxr = min(len(grid), 40)
    else:
        raise ValueError("method must be 'svd' or 'chebyshev'")
    cand = [rank] if rank else range(1, maxr + 1)
    for m in cand:
        fa, gb = build(m)
        err = float(np.max(np.abs(fa(grid) @ gb(grid).T - A) / A))
        if rank or err < tol:
            return fa, gb, m, err
    return fa, gb, m, err


class SeparableKernel:
    """terms: list of (kernel, a_expr, b_expr); a_expr (target weight) and
    b_expr (source weight) are sympy expressions in the symbols `params`."""

    def __init__(self, params, terms):
        self.params = tuple(params)
        self.terms = list(terms)
        keys = {}
        self.groups = []  # (kernel, [term index, ...])
        for t, (k, a, b) in enumerate(self.terms):
            key = (k.dim, str(k.expr))
            if key not in keys:
                keys[key] = len(self.groups)
                self.groups.append((k, []))
            self.groups[keys[key]][1].append(t)
        # weights may be sympy expressions in `params` or plain callables of the parameter columns
        self._a = [a if callable(a) else sp.lambdify(self.params, a, "numpy") for _, a, _ in self.terms]
        self._b = [b if callable(b) else sp.lambdify(self.params, b, "numpy") for _, _, b in self.terms]

    @property
    def nchannels(self):
        return len(self.terms)

    # -- constructors for Lennard-Jones --------------------------------------
    @classmethod
    def lennard_jones(cls, mixing="geometric", dim=3):
        """4 eps_ij ((sigma_ij/r)^12 - (sigma_ij/r)^6) with per-particle (eps, sigma)."""
        eps, sig = sp.symbols("eps sigma", positive=True)
        terms = []
        for n, sign in ((12, 1), (6, -1)):
            G = power_law(n, dim)
            if mixing == "geometric":
                terms.append((G, sign * 4 * sp.sqrt(eps) * sig ** (n // 2), sp.sqrt(eps) * sig ** (n // 2)))
            elif mixing in ("lorentz-berthelot", "lb"):
                for k in range(n + 1):
                    terms.append((G, sign * 4 * sp.sqrt(eps) * sp.Rational(math.comb(n, k), 2**n) * sig**k,
                                  sp.sqrt(eps) * sig ** (n - k)))
            else:
                raise ValueError("mixing must be 'geometric' or 'lorentz-berthelot'")
        return cls((eps, sig), terms)

    @classmethod
    def lennard_jones_lowrank(cls, sigma_sample=None, sigma_range=None, tol=1e-4, ranks=None,
                              method="svd", dim=3):
        """Lorentz-Berthelot Lennard-Jones with the sigma channels compressed to low rank.

        The exact rank of ((sigma_i + sigma_j)/2)^n is n+1, but sigma is bounded, so the pair matrix is
        numerically low rank. `method="svd"` builds the optimal basis from the sampled (or unique) sigma
        values by a Nystrom extension (weights are then evaluated for any sigma); `"chebyshev"` interpolates
        in sigma_j at Chebyshev nodes on `sigma_range` and needs no data. The rank per exponent is the
        smallest meeting the maximum entrywise relative error `tol` on the sigma grid, or `ranks=(r12, r6)`.
        epsilon_ij = sqrt(eps_i eps_j) stays exactly rank 1 and multiplies every channel.
        """
        eps, sig = sp.symbols("eps sigma", positive=True)
        terms, info = [], {}
        for n, sign, r in ((12, 1, None if ranks is None else ranks[0]), (6, -1, None if ranks is None else ranks[1])):
            fa, gb, rank, err = low_rank_sigma_factors(n, sigma_sample, sigma_range, tol, r, method)
            info[n] = (rank, err)
            G = power_law(n, dim)
            for k in range(rank):
                terms.append((G, (lambda e, s_, k=k, fa=fa, sign=sign: sign * 4 * np.sqrt(e) * fa(s_)[:, k]),
                              (lambda e, s_, k=k, gb=gb: np.sqrt(e) * gb(s_)[:, k])))
        sk = cls((eps, sig), terms)
        sk.compression = info
        return sk

    # -- weights -----------------------------------------------------------
    def target_weights(self, group, theta):
        """a_c for the channels of `group`, theta an (n, nparams) array -> (n, C)."""
        th = np.atleast_2d(theta)
        return np.stack([np.broadcast_to(self._a[t](*th.T), len(th)) for t in self.groups[group][1]], axis=1)

    def source_weights(self, group, theta):
        th = np.atleast_2d(theta)
        return np.stack([np.broadcast_to(self._b[t](*th.T), len(th)) for t in self.groups[group][1]], axis=1)

    # -- reference ---------------------------------------------------------
    def direct(self, targets, theta_t, sources, theta_s):
        """phi and F = -grad phi at targets by direct summation (numpy)."""
        d = targets.shape[1]
        phi = np.zeros(len(targets))
        F = np.zeros((len(targets), d))
        for t, (k, _, _) in enumerate(self.terms):
            G, grads = k.lambdify()
            a = np.broadcast_to(self._a[t](*np.atleast_2d(theta_t).T), len(targets))
            b = np.broadcast_to(self._b[t](*np.atleast_2d(theta_s).T), len(sources))
            for i, x in enumerate(targets):
                dd = x[None, :] - sources
                cols = [dd[:, j] for j in range(d)]
                phi[i] += a[i] * np.sum(b * G(*cols))
                for j in range(d):
                    F[i, j] -= a[i] * np.sum(b * grads[j](*cols))
        return phi, F

    # -- code generation ---------------------------------------------------
    def generate(self, name, order, orders=None, **kw):
        """One multi-channel operator set per group: {name}_g{g}.h/.c.

        Returns {group: report}. Channel order inside a group is that of
        `self.groups[g][1]`. Every group defines the same operator names, so each gets its own
        namespace: a function prefix `{name}_g{g}_` in C and C++ (all groups then link into one
        program) and, in Fortran, its module.
        """
        fortran = kw.get("language") == "fortran"
        return {g: generate_kernel_code(k, order, f"{name}_g{g}", orders=orders, nchannels=len(ts),
                                        prefix="" if fortran else f"{name}_g{g}_", **kw)
                for g, (k, ts) in enumerate(self.groups)}
