"""Tensor-valued kernels, described by an array of sympy expressions and an
einsum-style contraction with the sources.

    G = np.eye(3) / r + np.einsum("i,j->ij", X, X) / r**3       # Stokeslet
    k = Kernel(G, contract="ij,j->i")                            # -> VectorKernel

`contract` is "kernel,source->output": kernel axes that appear in the output are
target-side components, the others are summed against the source tensor (flattened
C-order over the source letters). It defaults to "ij,j->i" for a matrix kernel.

A VectorKernel is lowered to *groups*, each a scalar `Kernel` plus the sources it
takes, which the existing scalar machinery generates (`generate_kernel_code`):

* split_channels=True (and the fallback): entry E[o, s] is its own scalar channel
  fed with source component f_s and added to output o. Always correct; identical
  entries (G_ij = G_ji, equal diagonals) share one derivative tensor.
* a registered decomposition, tried by default when the constructor is given one
  (`xx_kernel`: a delta_ij / r + b x_i x_j / r^3 covers Stokes and Kelvin): the
  entries are differential operators on one scalar potential, so all outputs share
  one derivative tensor and one multipole array per component.

Generated operators return F = [phi, -dphi/dx_i] per channel, channel-major;
`output` folds that into the tensor output (and its gradient).
"""
from dataclasses import dataclass, field

import numpy as np
import sympy as sp

from .kernel import Kernel, generate_kernel_code, shell_indices

__all__ = ["VectorKernel", "Group", "xx_kernel", "stokeslet", "kelvin"]

_BIH3 = {(4, 0, 0): 1, (0, 4, 0): 1, (0, 0, 4): 1, (2, 2, 0): 2, (2, 0, 2): 2, (0, 2, 2): 2}


@dataclass
class Group:
    """One generated operator set: `kernel` with `nchannels` channels of order-`source_order` sources.

    out_comp[c] is the flat output component channel c adds into. weights(f) maps source
    tensors f (n, n_src) to the channel-major source coefficients (n, nchannels * shell).
    """

    kernel: Kernel
    source_order: int
    out_comp: list
    weights: callable = field(repr=False)

    @property
    def nchannels(self):
        return len(self.out_comp)


def _parse_contract(spec, ndim):
    if spec is None:
        if ndim != 2:
            raise ValueError("give contract='kernel,source->output' for kernels that are not matrices")
        spec = "ij,j->i"
    lhs, out = spec.replace(" ", "").split("->")
    ker, src = lhs.split(",")
    if len(ker) != ndim or len(set(ker)) != ndim:
        raise ValueError(f"contract {spec!r}: kernel needs {ndim} distinct indices")
    if set(out) - set(ker) or set(src) != set(ker) - set(out) or len(set(src)) != len(src):
        raise ValueError(f"contract {spec!r}: the source letters must be exactly the summed kernel letters")
    return ker, src, out


class VectorKernel:
    def __init__(self, array, variables=None, params=None, contract=None, split_channels=None,
                 decomposition=None, **kernel_kw):
        arr = np.empty(np.shape(array), dtype=object)
        flat = np.asarray(array, dtype=object)
        for idx in np.ndindex(arr.shape):
            arr[idx] = sp.sympify(flat[idx])
        self.array = arr
        if variables is None:
            names = {s.name: s for e in arr.flat for s in e.free_symbols}
            variables = tuple(names[n] for n in ("x", "y", "z") if n in names)
        self.variables = tuple(variables)
        self.dim = len(self.variables)
        self.params = params
        self.kernel_kw = kernel_kw
        ker, src, out = _parse_contract(contract, arr.ndim)
        self.contract = f"{ker},{src}->{out}"
        self.out_axes = [ker.index(c) for c in out]
        self.src_axes = [ker.index(c) for c in src]
        self.out_shape = tuple(arr.shape[a] for a in self.out_axes)
        self.src_shape = tuple(arr.shape[a] for a in self.src_axes)
        self.n_out = int(np.prod(self.out_shape))
        self.n_src = int(np.prod(self.src_shape))
        # E[o, s]: the kernel as a (n_out, n_src) matrix of scalar entries
        self.E = np.transpose(arr, self.out_axes + self.src_axes).reshape(self.n_out, self.n_src)
        if split_channels is None:
            split_channels = decomposition is None
        self.split_channels = bool(split_channels)
        if decomposition is not None and not self.split_channels:
            self.groups = decomposition
            self._check_groups()
        else:
            self.groups = self._split_groups()

    # -- lowering ---------------------------------------------------------
    def _split_groups(self):
        keys, entries = {}, []
        for o in range(self.n_out):
            for s in range(self.n_src):
                e = self.E[o, s]
                if e == 0:
                    continue
                key = sp.srepr(sp.simplify(e))
                keys.setdefault(key, []).append((o, s))
                if len(keys[key]) == 1:
                    entries.append((key, e))
        groups = []
        for key, e in entries:
            pairs = keys[key]
            k = Kernel(e, self.variables, params=self.params, **self.kernel_kw)
            cols = [s for _, s in pairs]
            groups.append(Group(k, 0, [o for o, _ in pairs],
                                lambda f, cols=cols: np.asarray(f, float).reshape(len(f), -1)[:, cols]))
        return groups

    def _check_groups(self):
        """Compare the decomposition with the array at a random point and force."""
        rng = np.random.default_rng(0)
        pt = rng.uniform(0.5, 1.5, self.dim) * rng.choice((-1, 1), self.dim)
        f = rng.standard_normal((1, self.n_src))
        ref = self._eval_entries(pt[None, :]) @ f[0]
        got = np.zeros(self.n_out)
        for g in self.groups:
            S = g.weights(f).reshape(g.nchannels, -1)
            sh = shell_indices(self.dim, g.source_order)
            for c, o in enumerate(g.out_comp):
                for m, im in sh.items():
                    d = g.kernel.expr
                    for v, e in zip(g.kernel.variables, m):
                        if e:
                            d = sp.diff(d, v, e)
                    got[o] += S[c, im] * float(d.subs(dict(zip(g.kernel.variables, pt))))
        if not np.allclose(got, ref, rtol=1e-9, atol=1e-12):
            raise ValueError(f"the decomposition does not reproduce the kernel (got {got}, expected {ref})")

    # -- numerics ----------------------------------------------------------
    def _eval_entries(self, d):
        """E at displacements d (n, dim): array (n, n_out, n_src)."""
        out = np.zeros((len(d), self.n_out, self.n_src))
        cols = [d[:, j] for j in range(self.dim)]
        for o in range(self.n_out):
            for s in range(self.n_src):
                e = self.E[o, s]
                if self.params:
                    e = e.subs({sp.Symbol(k) if isinstance(k, str) else k: v for k, v in self.params.items()})
                if e != 0:
                    out[:, o, s] = sp.lambdify(self.variables, e, "numpy")(*cols)
        return out

    def direct(self, targets, sources, f):
        """Output at each target by direct summation; f is (n_sources, n_src)."""
        t, y = np.asarray(targets, float), np.asarray(sources, float)
        f = np.asarray(f, float).reshape(len(y), -1)
        res = np.zeros((len(t), self.n_out))
        for i, ti in enumerate(t):
            res[i] = np.einsum("jos,js->o", self._eval_entries(ti - y), f)
        return res.reshape((len(t),) + self.out_shape)

    # -- operators ------------------------------------------------------------
    def generate(self, name, order, **kw):
        """{name}_g{g}.h/.c for each group; returns {group: report}."""
        # Every group defines the same operator names, so each gets its own namespace: a function
        # prefix in C and C++ (all groups then link into one program) and, in Fortran, its module.
        fortran = kw.get("language") == "fortran"
        return {g: generate_kernel_code(grp.kernel, order, f"{name}_g{g}", source_order=grp.source_order,
                                        nchannels=grp.nchannels, prefix="" if fortran else f"{name}_g{g}_", **kw)
                for g, grp in enumerate(self.groups)}

    def source_weights(self, g, f):
        return self.groups[g].weights(np.asarray(f, float).reshape(len(f), -1))

    def output(self, g, F, out, grad=None):
        """Add group g's generated-operator output F (one target) into out (n_out,) and grad (n_out, dim)."""
        grp = self.groups[g]
        Fc = np.asarray(F).reshape(grp.nchannels, 1 + self.dim)
        for c, o in enumerate(grp.out_comp):
            out[o] += Fc[c, 0]
            if grad is not None:
                grad[o] -= Fc[c, 1:]


def xx_kernel(a=1.0, b=1.0, **kw):
    """G_ij = a delta_ij / r + b x_i x_j / r^3 in 3D (Stokes: a = b = 1; Kelvin: a = 3 - 4 nu, b = 1).

    Uses x_i x_j / r^3 = delta_ij / r - d_i d_j r and 1/r = (Laplace r) / 2, so
    G f = [c f_i Laplace - b sum_j f_j d_i d_j] r with c = (a + b) / 2: one biharmonic
    kernel r, order-2 sources, three channels sharing one derivative tensor.
    """
    x, y, z = X = sp.symbols("x y z")
    r = sp.sqrt(x**2 + y**2 + z**2)
    a, b = sp.nsimplify(a), sp.nsimplify(b)
    G = np.array([[(a if i == j else 0) / r + b * X[i] * X[j] / r**3 for j in range(3)] for i in range(3)],
                 dtype=object)
    shell = shell_indices(3, 2)
    c = (a + b) / 2

    def weights(f):
        f = np.asarray(f, float)
        S = np.zeros((len(f), 3, len(shell)))
        for m, im in shell.items():
            if max(m) == 2:
                k = m.index(2)
                for i in range(3):
                    S[:, i, im] = float(c) * f[:, i] - (float(b) * f[:, i] if k == i else 0.0)
            else:
                p, q = [u for u, e in enumerate(m) if e]
                S[:, p, im] = -float(b) * f[:, q]
                S[:, q, im] = -float(b) * f[:, p]
        return S.reshape(len(f), -1)

    kernel = Kernel(r, X, pde=_BIH3)
    return VectorKernel(G, X, contract="ij,j->i", decomposition=[Group(kernel, 2, [0, 1, 2], weights)], **kw)


def stokeslet(**kw):
    """delta_ij / r + x_i x_j / r^3 (viscosity 1 / (8 pi))."""
    return xx_kernel(1.0, 1.0, **kw)


def kelvin(nu=0.25, **kw):
    """(3 - 4 nu) delta_ij / r + x_i x_j / r^3: Kelvin's point force in an elastic solid,
    up to 1 / (16 pi mu (1 - nu))."""
    return xx_kernel(3 - 4 * sp.nsimplify(nu), 1.0, **kw)
