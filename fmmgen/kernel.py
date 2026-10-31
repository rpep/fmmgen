"""
Kernel-generic Cartesian FMM operators.

The original fmmgen path is hard-wired to the Laplace kernel 1/r in 3D. This
module lifts that: describe a translation-invariant kernel G(r) once, with
`Kernel`, and `generate_kernel_code` emits S2M/M2M/M2L/L2L/L2P/M2P/P2P
(potential, optionally field) through the existing printers, in C, C++ or
Fortran 90.

Conventions (identical to the Laplace path so the two are interchangeable):

    phi(t) = sum_s q_s G(t - s)
    M_n    = sum_s q_s (c - s)^n / n!             multipole about centre c
    L_n    = sum_m M_m D_{n+m}(R),  R = t_c - c   local about centre t_c
    phi(t) = sum_n L_n (t - t_c)^n / n!
    D_n(R) = d^n G / dr^n at r = R

M2M and L2L are Taylor shifts and do not depend on the kernel. The kernel
enters only through the derivative tensor D_n, and that is the part with
several strategies (see `Kernel` and `derivative_layers`):

    "recurrence" radial kernels, the default. G = g(R) = f(s), s = R^2. With
                E[n][k] = d^n f^(k)(s) one has d_i f^(k) = 2 x_i f^(k+1), so by
                Leibniz  E[n+e_i][k] = 2 x_i E[n][k+1] + 2 n_i E[n-e_i][k+1]
                and D_n = E[n][0]. Only |n| + k <= p is needed: ~3 flops per
                entry of a C(p+d+1, d+1) table (O(p^(d+1)); 495 entries at
                d=3, p=8) after p univariate derivatives f^(k) -- no
                multivariate symbolic differentiation, no expression swell.
    "faa"       the same tensor in closed form: n! [t^n] f(s + 2R.t + t.t)
                expanded by multivariate Faa di Bruno, i.e. D_n = sum_k
                f^(k)/k! * n! [t^n] (2R.t + t.t)^k. Kept for comparison: each
                D_n is an independent polynomial, so the generated code is
                larger than the recurrence's.
    "symbolic"  sp.diff of the full expression, reusing D_{n-e_i}. Works for any
                (including non-radial) kernel but the expression swell is real:
                see research/kernels/FINDINGS.md.
    pde         optional. A constant-coefficient linear PDE the kernel obeys
                away from the origin (Laplace, Helmholtz, Yukawa, ...) turns the
                derivatives with n_last >= pivot degree into recurrences on
                lower ones, so only a minority are evaluated by the radial scheme. Pass
                `pde={(2,0,0): 1, (0,2,0): 1, (0,0,2): 1}`.

Anything else (e.g. the PDE-compression work) can plug in through the `derivs`
hook: a callable `derivs(kernel, pmax, index_map) -> layers`, where layers is a
list of lists of `(n, expr)`. `expr` may use the kernel's variables, the
symbol `R` and `Dref(m)` for the derivative
with multi-index m defined in an EARLIER layer.
"""
import itertools
import math
import time

import sympy as sp
from sympy import ring, ZZ
from sympy.matrices.expressions.matexpr import MatrixElement

from .printers import FunctionPrinter

__all__ = ["Kernel", "generate_kernel_code", "multi_indices", "Dref", "derivative_layers",
           "numeric_derivative_tensor", "power_law", "lennard_jones", "morse", "buckingham", "ewald_real",
           "gaussian_screened_coulomb", "softened_power", "yukawa", "helmholtz_cos", "helmholtz_sin"]

R = sp.Symbol("R")  # plain symbol: the printer defines it when it sees it
_Rpos = sp.Symbol("R", positive=True)


# ----------------------------------------------------------------------------
# Multi-index helpers (any dimension)
# ----------------------------------------------------------------------------
def nterms(dim, p):
    """Number of multi-indices with |n| <= p in `dim` dimensions."""
    return math.comb(p + dim, dim) if p >= 0 else 0


def multi_indices(dim, p):
    """{n: index} for |n| <= p, ordered like utils.generate_mappings (grevlex).

    For dim=3 this is exactly the existing Laplace array layout, so arrays are
    interchangeable with the original operators.
    """
    ns = [n for n in itertools.product(range(p + 1), repeat=dim) if sum(n) <= p]
    ns.sort(key=lambda n: (sum(n), tuple(-e for e in n)))
    return {n: i for i, n in enumerate(ns)}


def _fact(n):
    out = 1
    for e in n:
        out *= math.factorial(e)
    return out


def _add(a, b):
    return tuple(i + j for i, j in zip(a, b))


def _sub(a, b):
    return tuple(i - j for i, j in zip(a, b))


def _monomial(vars_, n):
    out = sp.Integer(1)
    for v, e in zip(vars_, n):
        out *= v**e
    return out


def Dref(m):
    """Placeholder for the derivative D_m defined in an earlier layer."""
    return sp.Symbol("Dref_" + "_".join(str(e) for e in m))


def _as_number(v):
    if isinstance(v, float):
        return sp.Rational(repr(v))
    return sp.sympify(v)


FK = "fk"  # name of the univariate-derivative array in generated code


# ----------------------------------------------------------------------------
# Kernel description
# ----------------------------------------------------------------------------
class Kernel:
    """A translation-invariant kernel G(r), r = target - source.

    Parameters
    ----------
    expr : sympy expression in `variables` (and `params`).
    variables : tuple of symbols; its length is the dimension (1, 2 or 3).
        Defaults to the free symbols x, y, z that appear, in that order.
    params : {symbol or name: number}. Numerical values are baked into the
        generated code (regenerate to change them).
    pde : {multi-index: coefficient} describing  sum c_a d^a G = 0  (r != 0).
        The highest degree in the last variable must occur in exactly one term,
        with all other terms of that degree lower, and no term may have higher
        total order than the pivot. It is verified numerically at construction.
    derivs : callable hook, see the module docstring.
    method : "auto" (default), "recurrence", "faa" or "symbolic".
    radial_derivs : optional callable kmax -> [f_0..f_kmax] (sympy in the radial
        symbol), the univariate derivatives f^(k)(s) of the profile in s = R^2,
        for kernels whose naive derivatives cancel catastrophically (see
        `_stable_erf_derivs`). Default: differentiate the profile symbolically.

    Radial kernels (G depends on |r| only) are detected numerically; their
    radial profile g(R) is then available as `.profile` and enables the
    "taylor" method and a cheap closed-form P2P field.
    """

    def __new__(cls, expr=None, *args, **kw):
        # an array-valued expression (numpy object array, nested list/tuple, sympy Array) is a tensor kernel
        if cls is Kernel and (isinstance(expr, (list, tuple, sp.Array)) or getattr(expr, "ndim", 0) > 0):
            from .vector import VectorKernel

            return VectorKernel(expr, *args, **kw)
        return super().__new__(cls)

    def __init__(self, expr, variables=None, params=None, pde=None, derivs=None, method="auto", name=None,
                 radial_derivs=None):
        expr = sp.sympify(expr)
        if variables is None:
            names = {s.name: s for s in expr.free_symbols}
            variables = tuple(names[n] for n in ("x", "y", "z") if n in names)
        self.variables = tuple(variables)
        self.dim = len(self.variables)
        if not 1 <= self.dim <= 3:
            raise ValueError("dimension must be 1, 2 or 3")
        subs = {}
        for k, v in (params or {}).items():
            k = k if isinstance(k, sp.Symbol) else sp.Symbol(k)
            subs[k] = _as_number(v)
        self.params = subs
        expr = expr.subs(subs)
        leftover = expr.free_symbols - set(self.variables)
        if leftover:
            raise ValueError(f"unbound symbols {sorted(map(str, leftover))}: give values via params=")
        self.expr = expr
        self.name = name or str(expr)
        self.pde = dict(pde) if pde else None
        self.derivs = derivs
        self._radial_derivs_fn = radial_derivs
        if method not in ("auto", "recurrence", "faa", "symbolic"):
            raise ValueError("method must be auto, recurrence, faa or symbolic")
        self.profile = self._radial_profile()
        if method == "auto":
            method = "recurrence" if self.profile is not None else "symbolic"
        if method in ("recurrence", "faa") and self.profile is None:
            raise ValueError("method='recurrence'/'faa' need a radial kernel")
        self.method = method
        if self.pde:
            self._check_pde()

    @classmethod
    def radial(cls, g, r, dim=3, variables=None, **kw):
        """Kernel from a profile g(r) written in a symbol `r`."""
        if variables is None:
            variables = sp.symbols("x y z")[:dim]
        r2 = sum(v**2 for v in variables)
        return cls(sp.sympify(g).subs(r, sp.sqrt(r2)), variables, **kw)

    # -- radial detection ------------------------------------------------
    def _radial_profile(self):
        v = self.variables
        g = self.expr.subs({v[0]: _Rpos, **{w: 0 for w in v[1:]}})
        if g.has(sp.zoo, sp.nan, sp.oo):
            return None
        import random

        rng = random.Random(1)
        for _ in range(4):
            pt = [rng.uniform(0.4, 1.7) * rng.choice((-1, 1)) for _ in v]
            r = math.sqrt(sum(c * c for c in pt))
            try:
                a = complex(self.expr.subs(dict(zip(v, pt))).evalf())
                b = complex(g.subs(_Rpos, r).evalf())
            except TypeError:
                return None
            if abs(a - b) > 1e-9 * (1 + abs(a)):
                return None
        return g

    def _check_pde(self):
        v = self.variables
        pt = dict(zip(v, [0.7, -0.4, 0.9][: self.dim]))
        total = 0
        scale = 0
        for a, c in self.pde.items():
            if len(a) != self.dim:
                raise ValueError("pde multi-index length must equal the dimension")
            d = self.expr
            for var, e in zip(v, a):
                if e:
                    d = sp.diff(d, var, e)
            term = complex(_as_number(c).subs(self.params) * d.subs(pt).evalf())
            total += term
            scale += abs(term)
        if abs(total) > 1e-8 * (1 + scale):
            raise ValueError(f"the kernel does not satisfy the supplied pde (residual {abs(total):.2e})")

    def radial_derivs(self, kmax):
        """f^(k)(s), k = 0..kmax, as expressions in the radial symbol."""
        if self._radial_derivs_fn is not None:
            return self._radial_derivs_fn(kmax)
        return _radial_derivs(self.profile, kmax)

    # -- numerics for tests ----------------------------------------------
    def lambdify_expr(self, expr):
        """numpy function of the kernel variables for `expr` (erf/erfc vectorised from math)."""
        import numpy as np

        mods = [{"erf": np.vectorize(math.erf), "erfc": np.vectorize(math.erfc)}, "numpy"]
        return sp.lambdify(self.variables, expr, mods)

    def lambdify(self):
        """numpy function G(*coords) and its gradient, for direct summation."""
        G = self.lambdify_expr(self.expr)
        grads = [self.lambdify_expr(sp.diff(self.expr, w)) for w in self.variables]
        return G, grads

    def __repr__(self):
        return f"Kernel({self.expr}, dim={self.dim}, method={self.method!r}, pde={bool(self.pde)})"


# ----------------------------------------------------------------------------
# Derivative tensor
# ----------------------------------------------------------------------------
def _radial_derivs(g, kmax):
    """f_k(R) = f^(k)(s) at s = R^2, k=0..kmax, for the profile g(R) = f(sqrt(s)).

    Differentiated in s itself (a single-variable problem). Doing it through
    ((1/(2R)) d/dR)^k instead leaves unsimplified R/R factors that swell badly
    for kernels like (R^2 + eps^2)^(-1/2); measured in FINDINGS.md.
    """
    s = sp.Symbol("s", positive=True)
    f = g.subs(_Rpos, sp.sqrt(s))
    fs = [f]
    for _ in range(kmax):
        fs.append(sp.diff(fs[-1], s))
    return [sp.powsimp(sp.expand(e.subs(s, _Rpos**2))) for e in fs]


def _stable_erf_derivs(alpha, kmax, extra=60):
    """f^(k)(s) for f(s) = erf(alpha sqrt(s)) / sqrt(s), free of cancellation.

    f(s) = (2 alpha / sqrt(pi)) int_0^1 exp(-alpha^2 s t^2) dt, so
    f^(k)(s) = (2 alpha / sqrt(pi)) (-alpha^2)^k I_k(c), c = alpha^2 s, with the
    positive integrals I_k = int_0^1 t^(2k) exp(-c t^2) dt. The naive closed form
    is a sum of terms growing like k! that cancel to a small result, losing ~k
    digits (all of float64 by k = 17). I_k obeys
        I_k = ((2k-1) I_(k-1) - e^-c) / (2c)    stable upward for c > kmax,
        I_(k-1) = (2c I_k + e^-c) / (2k-1)      all terms positive, so stable
                                                downward from k = kmax + extra.
    The result is a Piecewise on c, accurate to ~1e-15 for c = 0.001..6400, k <= 24.
    """
    # All coefficients are 60-digit Floats, never Rationals: with exact rationals sympy folds
    # the chain into Horner polynomials whose integer coefficients (~1e200) overflow C literals.
    def F(v):
        return sp.Float(v, 60)

    alpha = F(_as_number(alpha).evalf(60))
    c = alpha**2 * _Rpos**2
    e = sp.exp(-c)
    up = [sp.sqrt(sp.pi) * sp.erf(alpha * _Rpos) / (2 * alpha * _Rpos)]
    for k in range(1, kmax + 1):
        up.append((F(2 * k - 1) * up[-1] - e) / (2 * c))
    n = kmax + extra
    dn = {n: e / F(2 * n + 1)}
    for k in range(n, 0, -1):
        dn[k - 1] = (2 * c * dn[k] + e) / F(2 * k - 1)
    def pref(k):
        return (2 * alpha / sp.sqrt(sp.pi) * (-(alpha**2)) ** k).evalf(60)

    return [sp.Piecewise((pref(k) * up[k], c > kmax), (pref(k) * dn[k], True)) for k in range(kmax + 1)]


def _delta_power_coeffs(dim, pmax):
    """C[k][n] = [t^n] (2 x.t + t.t)^k as ring elements in x (n: |n|<=pmax).

    Returns (ring, C). C[k] is a dict n -> PolyElement. The univariate
    derivative f_k enters D_n for ceil(|n|/2) <= k <= |n|.
    """
    names = ",".join(f"X{i}" for i in range(dim))
    Rg, *X = ring(names, ZZ)
    delta = {}
    for i in range(dim):
        e1 = tuple(1 if j == i else 0 for j in range(dim))
        e2 = tuple(2 if j == i else 0 for j in range(dim))
        delta[e1] = delta.get(e1, Rg.zero) + 2 * X[i]
        delta[e2] = delta.get(e2, Rg.zero) + Rg.one
    zero_idx = tuple([0] * dim)
    C = [{zero_idx: Rg.one}]
    for k in range(1, pmax + 1):
        new = {}
        for n1, c1 in C[-1].items():
            d1 = sum(n1)
            for n2, c2 in delta.items():
                if d1 + sum(n2) <= pmax:
                    n = _add(n1, n2)
                    new[n] = new.get(n, Rg.zero) + c1 * c2
        C.append(new)
    return Rg, C


def _taylor_entries(kernel, pmax, ns):
    """Expressions for D_n, n in ns, in terms of x_i and fk[k]."""
    dim, vars_ = kernel.dim, kernel.variables
    Rg, C = _delta_power_coeffs(dim, pmax)
    fk = sp.MatrixSymbol(FK, pmax + 1, 1)
    out = []
    for n in ns:
        deg = sum(n)
        acc = sp.Integer(0)
        for k in range(math.ceil(deg / 2), deg + 1):
            c = C[k].get(n)
            if c is None or c == 0:
                continue
            poly = c.as_expr(*vars_)
            acc += sp.Rational(_fact(n), math.factorial(k)) * poly * fk[k]
        out.append((n, acc))
    return out


def _to_plain_R(e):
    return e.xreplace({_Rpos: R})


def _symbolic_layer(kernel, pmax):
    """Per-entry sp.diff, reusing D_{n-e_i}; one layer, expressed with R."""
    vars_ = kernel.variables
    idx = multi_indices(kernel.dim, pmax)
    r2 = sum(v**2 for v in vars_)
    D = {}
    entries = []
    for n in idx:
        if sum(n) == 0:
            e = kernel.expr
        else:
            i = next(j for j, ex in enumerate(n) if ex)
            prev = _sub(n, tuple(1 if j == i else 0 for j in range(kernel.dim)))
            e = sp.diff(D[prev], vars_[i])
        D[n] = e
    for n in idx:
        e = D[n]
        if kernel.profile is not None:
            e = e.replace(lambda t: t.is_Pow and t.base == r2, lambda t: R ** (2 * t.exp))
        entries.append((n, e))
    return entries


def _eref(n, k):
    return sp.Symbol("Dref_E_" + "_".join(str(e) for e in n) + f"_{k}")


def _first_nonzero(n):
    return next(i for i, e in enumerate(n) if e)


def _recurrence_layers(kernel, pmax, targets):
    """Radial-function recurrence in s = r^2 for D_n, n in `targets`.

    Write E[n][k] = d^n f^(k)(s) (f^(k) the k-th s-derivative of the profile
    f(s) = g(sqrt(s))). Since d_i f^(k) = 2 x_i f^(k+1), Leibniz gives
        E[n+e_i][k] = 2 x_i E[n][k+1] + 2 n_i E[n-e_i][k+1],
    and D_n = E[n][0]. Only (n, k) with |n| + k <= p are needed, so the whole
    tensor costs ~3 flops per entry of a C(p+d+1, d+1)-sized table: O(p^(d+1)),
    against a p^d-sized output and the exponential swell of sp.diff. Layer j of
    the output holds the |n| = j entries; layer 0 holds the univariate
    f^(k)(R^2), the only place the kernel formula (and R) appears.
    Returns layers of (key, expr) with key = (\"E\", n, k).
    """
    dim, vars_ = kernel.dim, kernel.variables
    need = {}
    stack = [(n, 0) for n in targets]
    while stack:
        n, k = stack.pop()
        if k in need.setdefault(n, set()):
            continue
        need[n].add(k)
        if sum(n) > 0:
            i = _first_nonzero(n)
            m = tuple(e - (j == i) for j, e in enumerate(n))
            stack.append((m, k + 1))
            if m[i] > 0:
                stack.append((tuple(e - (j == i) for j, e in enumerate(m)), k + 1))
    kmax = max(need[tuple([0] * dim)])
    fs = kernel.radial_derivs(kmax)
    zero = tuple([0] * dim)
    layers = [[(("E", zero, k), _to_plain_R(fs[k])) for k in sorted(need[zero])]]
    for j in range(1, pmax + 1):
        layer = []
        for n in sorted((n for n in need if sum(n) == j), key=lambda n: multi_indices(dim, pmax)[n]):
            i = _first_nonzero(n)
            m = tuple(e - (t == i) for t, e in enumerate(n))
            mm = tuple(e - (t == i) for t, e in enumerate(m)) if m[i] > 0 else None
            for k in sorted(need[n]):
                e = 2 * vars_[i] * _eref(m, k + 1)
                if mm is not None:
                    e += 2 * m[i] * _eref(mm, k + 1)
                layer.append((("E", n, k), e))
        if layer:
            layers.append(layer)
    return layers


def _pde_split(kernel):
    pde = {a: _as_number(c).subs(kernel.params) for a, c in kernel.pde.items() if c != 0}
    zdeg = max(a[-1] for a in pde)
    piv = [a for a in pde if a[-1] == zdeg]
    rest = {a: c for a, c in pde.items() if a not in piv}
    if len(piv) != 1 or zdeg == 0 or any(a[-1] >= zdeg for a in rest) or any(
        sum(a) > sum(piv[0]) for a in rest
    ):
        raise NotImplementedError(
            "pde needs a unique pivot term of strictly highest degree in the last variable and no "
            "higher-order terms (Laplace, Helmholtz, Yukawa, biharmonic qualify); use the derivs= hook otherwise"
        )
    return piv[0], pde[piv[0]], rest, zdeg


def derivative_layers(kernel, pmax):
    """Layers of (key, expr) defining D_n for |n| <= pmax; see module docstring.

    key is a multi-index n (defines D_n) or ("E", n, k) (an auxiliary table
    entry; D_n is ("E", n, 0)). Returns (layers, names).
    """
    idx = multi_indices(kernel.dim, pmax)
    if kernel.derivs is not None:
        layers = kernel.derivs(kernel, pmax, idx)
    elif kernel.method == "symbolic":
        layers = [_symbolic_layer(kernel, pmax)]
    else:
        if kernel.method == "faa":
            mk = lambda ns: [_taylor_entries(kernel, pmax, ns)]  # noqa: E731
        else:
            mk = lambda ns: _recurrence_layers(kernel, pmax, ns)  # noqa: E731
        if not kernel.pde:
            layers = mk(list(idx))
        else:
            pivot, cpiv, rest, zdeg = _pde_split(kernel)
            base = [n for n in idx if n[-1] < zdeg]
            layers = mk(base)
            # D_n for base n is the layer's own entry; give it a D-name too
            layer_of = {n: 0 for n in base}
            pending = []
            for n in sorted((n for n in idx if n[-1] >= zdeg), key=lambda n: (n[-1], sum(n))):
                refs = {_add(_sub(n, pivot), a): c for a, c in rest.items()}
                expr = -sum((c / cpiv) * Dref(m) for m, c in refs.items())
                layer_of[n] = 1 + max(layer_of[m] for m in refs)
                pending.append((layer_of[n], n, expr))
            for L in range(1, max(lay for lay, _, _ in pending) + 1 if pending else 1):
                layers.append([(n, e) for lay, n, e in pending if lay == L])
    names = [f"D{j}" for j in range(len(layers))]
    return layers, names


class _DTensor:
    """Symbolic access D[n] to a derivative array built from layers."""

    def __init__(self, kernel, pmax):
        self.kernel = kernel
        self.pmax = pmax
        layers, names = derivative_layers(kernel, pmax)
        self.loc = {}
        self.arrays = []  # (name, Matrix)
        refs = {}
        uses_fk = False
        for lname, layer in zip(names, layers):
            sym = sp.MatrixSymbol(lname, len(layer), 1)
            exprs = []
            for i, (n, e) in enumerate(layer):
                if n and n[0] == "E":
                    refs[_eref(n[1], n[2])] = sym[i]
                    if n[2] == 0:
                        self.loc[n[1]] = sym[i]
                        refs[Dref(n[1])] = sym[i]
                else:
                    self.loc[n] = sym[i]
                    refs[Dref(n)] = sym[i]
            for n, e in layer:
                e = sp.sympify(e)
                e = e.xreplace(refs)
                uses_fk = uses_fk or any(isinstance(a, MatrixElement) and a.parent.name == FK
                                         for a in e.atoms(MatrixElement))
                exprs.append(e)
            self.arrays.append((lname, sp.Matrix(exprs)))
        if uses_fk:
            fs = kernel.radial_derivs(pmax)
            self.arrays.insert(0, (FK, sp.Matrix([_to_plain_R(f) for f in fs])))
        missing = [n for n in multi_indices(kernel.dim, pmax) if n not in self.loc]
        if missing:
            raise ValueError(f"derivative scheme did not define D_n for {missing[:3]}...")

    def __getitem__(self, n):
        return self.loc[n]


def numeric_derivative_tensor(kernel, pmax, point, mp_dps=None):
    """All D_n(point), |n| <= pmax, evaluated numerically by the radial recurrence.

    Same scheme as the generated code (see `_recurrence_layers`) but run at
    call time, so it is usable at orders where emitting and compiling straight
    line code is impractical (p = 16, 24, ...). Only the univariate radial
    derivatives f^(k) are symbolic (a one-variable problem); everything else is
    the O(p^(d+1)) table recurrence. With `mp_dps` set the arithmetic is done in
    mpmath at that many digits (to quantify cancellation); otherwise float64.
    Returns {n: value}.
    """
    if kernel.profile is None:
        raise ValueError("numeric_derivative_tensor needs a radial kernel")
    dim = kernel.dim
    fs = kernel.radial_derivs(pmax)
    if mp_dps:
        import mpmath

        mpmath.mp.dps = mp_dps
        x_ = [mpmath.mpf(c) for c in point]
        r = mpmath.sqrt(sum(c * c for c in x_))
        fk = [sp.lambdify(_Rpos, f, "mpmath")(r) for f in fs]
    else:
        x_ = [float(c) for c in point]
        r = math.sqrt(sum(c * c for c in x_))
        fk = [sp.lambdify(_Rpos, f, "math")(r) for f in fs]
    zero = tuple([0] * dim)
    E = {zero: fk}  # E[n][k], k = 0..pmax-|n|
    for j in range(1, pmax + 1):
        for n in (n for n in multi_indices(dim, j) if sum(n) == j):
            i = _first_nonzero(n)
            m = tuple(e - (t == i) for t, e in enumerate(n))
            Em = E[m]
            mm = tuple(e - (t == i) for t, e in enumerate(m)) if m[i] else None
            Emm = E[mm] if mm is not None else None
            c = 2 * x_[i]
            row = [c * Em[k + 1] for k in range(pmax - j + 1)]
            if Emm is not None:
                d = 2 * m[i]
                row = [row[k] + d * Emm[k + 1] for k in range(pmax - j + 1)]
            E[n] = row
    return {n: row[0] for n, row in E.items()}


# ----------------------------------------------------------------------------
# Operator expressions (source order s, nch channels sharing one kernel)
#
# A source of order s sits at y with coefficients S_m, |m| = s, and produces
#     phi(x) = sum_{|m|=s} S_m  d^m G (x - y)
# (derivatives with respect to the argument; a physical dipole mu has
# S = -mu because (mu . grad_y) G(x - y) = -(mu . grad) G). Its multipole about
# c is M_n = sum_m S_m (c - y)^(n-m) / (n-m)!  for s <= |n| <= p, and the local
# expansion has order p - s (|n| <= p - s), as in the Laplace source_order path.
# `nch` channels are stacked channel-major in S, M, L and the outputs; they share
# the derivative tensor, which is what makes mixing-rule kernels cheap.
# ----------------------------------------------------------------------------
def source_indices(dim, p, s):
    """{n: index} for the multipole array: s <= |n| <= p, legacy ordering."""
    ns = [n for n in multi_indices(dim, p) if sum(n) >= s]
    return {n: i for i, n in enumerate(ns)}


def shell_indices(dim, s):
    return source_indices(dim, s, s)


def _ops_S2M(dim, vars_, p, s=0, nch=1):
    J, sh = source_indices(dim, p, s), shell_indices(dim, s)
    S = sp.MatrixSymbol("S", len(sh) * nch, 1)
    out = []
    for c in range(nch):
        for n in J:
            acc = sp.Integer(0)
            for m, im in sh.items():
                k = _sub(n, m)
                if min(k) >= 0:
                    acc += S[c * len(sh) + im] * _monomial(vars_, k) / _fact(k)
            out.append(acc)
    return out


def _ops_M2M(dim, vars_, p, s=0, nch=1):
    J = source_indices(dim, p, s)
    nm = len(J)
    M = sp.MatrixSymbol("M", nm * nch, 1)
    out = []
    ks = multi_indices(dim, p - s)
    for c in range(nch):
        for n in J:
            acc = sp.Integer(0)
            for k in ks:
                nk = _sub(n, k)
                if nk in J:
                    acc += M[c * nm + J[nk]] * _monomial(vars_, k) / _fact(k)
            out.append(acc)
    return out


def _ops_L2L(dim, vars_, p, s=0, nch=1):
    pl = p - s
    idx = multi_indices(dim, pl)
    nl = len(idx)
    L = sp.MatrixSymbol("L", nl * nch, 1)
    out = []
    for c in range(nch):
        for n in idx:
            acc = sp.Integer(0)
            for k in multi_indices(dim, pl - sum(n)):
                acc += L[c * nl + idx[_add(n, k)]] * _monomial(vars_, k) / _fact(k)
            out.append(acc)
    return out


def _ops_L2P(dim, vars_, p, potential, field, s=0, nch=1):
    idx = multi_indices(dim, p - s)
    nl = len(idx)
    L = sp.MatrixSymbol("L", nl * nch, 1)
    out = []
    for c in range(nch):
        phi = sum((L[c * nl + i] * _monomial(vars_, n) / _fact(n) for n, i in idx.items()), sp.Integer(0))
        if potential:
            out.append(phi)
        if field:
            out += [-sp.diff(phi, v) for v in vars_]
    return out


def _ops_M2L(dim, p, D, s=0, nch=1):
    J = source_indices(dim, p, s)
    nm = len(J)
    idx = multi_indices(dim, p - s)
    M = sp.MatrixSymbol("M", nm * nch, 1)
    out = []
    for c in range(nch):
        for n in idx:
            acc = [M[c * nm + j] * D[_add(n, m)] for m, j in J.items() if sum(n) + sum(m) <= p]
            out.append(sp.Add(*acc, evaluate=False) if acc else sp.Integer(0))
    return out


def _ops_M2P(dim, p, D, potential, field, s=0, nch=1):
    J = source_indices(dim, p, s)
    nm = len(J)
    M = sp.MatrixSymbol("M", nm * nch, 1)
    out = []
    for c in range(nch):
        if potential:
            out.append(sp.Add(*[M[c * nm + j] * D[m] for m, j in J.items()], evaluate=False))
        if field:
            for d_ in range(dim):
                e = tuple(1 if k == d_ else 0 for k in range(dim))
                out.append(-sp.Add(*[M[c * nm + j] * D[_add(m, e)] for m, j in J.items()], evaluate=False))
    return out


def _flat_D(kernel, pmax):
    """{n: expression in x_i, R} for the whole tensor, layers inlined.

    Only for small pmax (P2P with source order s needs s, or s+1 with field):
    P2P_batch emits a single array, so there is no room for internal arrays.
    """
    Dt = _DTensor(kernel, pmax)
    sub = {}
    for name, mat in Dt.arrays:
        sym = sp.MatrixSymbol(name, len(mat), 1)
        for i, e in enumerate(mat):
            sub[sym[i]] = sp.sympify(e).xreplace(sub)
    return {n: sub[loc] for n, loc in Dt.loc.items()}


def _ops_P2P(kernel, potential, field, s=0, nch=1):
    sh = shell_indices(kernel.dim, s)
    S = sp.MatrixSymbol("S", len(sh) * nch, 1)
    vars_ = kernel.variables
    out = []
    if s == 0 and kernel.profile is not None:
        g = _to_plain_R(kernel.profile)
        dg = _to_plain_R(sp.diff(kernel.profile, _Rpos))
        for c in range(nch):
            if potential:
                out.append(S[c] * g)
            if field:
                out += [-S[c] * dg * v / R for v in vars_]
        return out
    if s == 0:
        G = kernel.expr
        for c in range(nch):
            if potential:
                out.append(S[c] * G)
            if field:
                out += [-S[c] * sp.diff(G, v) for v in vars_]
        return out
    D = _flat_D(kernel, s + (1 if field else 0))
    for c in range(nch):
        if potential:
            out.append(sum((S[c * len(sh) + i] * D[m] for m, i in sh.items()), sp.Integer(0)))
        if field:
            for d_ in range(kernel.dim):
                e = tuple(1 if k == d_ else 0 for k in range(kernel.dim))
                out.append(-sum((S[c * len(sh) + i] * D[_add(m, e)] for m, i in sh.items()), sp.Integer(0)))
    return out


# ----------------------------------------------------------------------------
# Code generation
# ----------------------------------------------------------------------------


def _dedupe_R(code):
    """Drop repeated `R`/`Rinv` declarations from one C function body."""
    seen, out = set(), []
    for line in code.split("\n"):
        s = line.strip()
        if s.startswith(("double R =", "double Rinv =", "float R =", "float Rinv =")):
            if s in seen:
                continue
            seen.add(s)
        out.append(line)
    return "\n".join(out)


# ----------------------------------------------------------------------------
# Planar variants: sources and targets confined to the z = 0 plane (moments stay full 3D).
#
# Every displacement then has z = 0. That kills every z power in the Taylor shifts, so the multipole
# array can only have n_z <= source_order non-zero entries, and L2P at dz = 0 only reads local entries
# with n_z <= 1 (n_z <= 0 without the field). The planar operators are the full operators restricted to
# those entries, with z set to 0, so they are valid for any kernel.
# ----------------------------------------------------------------------------
def _planar_keep(dim, p, s, nch, field):
    """Global indices of the multipole and local entries kept by the planar variants (channel-major)."""
    J, idx = source_indices(dim, p, s), multi_indices(dim, p - s)
    zmax = 1 if field else 0
    keep_m = [c * len(J) + j for c in range(nch) for n, j in sorted(J.items(), key=lambda t: t[1]) if n[-1] <= s]
    keep_l = [c * len(idx) + i for c in range(nch) for n, i in sorted(idx.items(), key=lambda t: t[1])
              if n[-1] <= zmax]
    return keep_m, keep_l


def planar_sizes(dim, p, s, nch, field):
    """(multipole, local) array lengths of the planar variants at order p."""
    km, kl = _planar_keep(dim, p, s, nch, field)
    return len(km), len(kl)


def _planar_ops(ops, rows, inputs, z):
    """`ops[r] for r in rows` with z = 0 and the input arrays restricted to their kept entries.

    inputs: {array name: (full length, kept global indices)}; entries that are dropped become 0."""
    rep = {z: sp.Integer(0)}
    for name, (n, keep) in inputs.items():
        full, new = sp.MatrixSymbol(name, n, 1), sp.MatrixSymbol(name, len(keep), 1)
        pos = {g: k for k, g in enumerate(keep)}
        for g in range(n):
            rep[full[g]] = new[pos[g]] if g in pos else sp.Integer(0)
    return [sp.sympify(ops[r]).xreplace(rep) for r in rows]


def _planar_arrays(arrays, ops, z):
    """Set z = 0 in the internal derivative arrays, drop the entries that vanish or are never used, and
    renumber the rest. Returns (arrays, ops)."""
    zero = sp.Integer(0)
    rep = {z: zero}
    sym = {name: sp.MatrixSymbol(name, len(mat), 1) for name, mat in arrays}
    exprs = {}
    for name, mat in arrays:
        exprs[name] = []
        for i, e in enumerate(mat):
            e = sp.sympify(e).xreplace(rep)
            exprs[name].append(e)
            if e == 0:
                rep[sym[name][i]] = zero  # every later reference to a vanishing entry vanishes too
    ops = [e.xreplace(rep) for e in ops]

    def used(e):
        return {a for a in e.atoms(MatrixElement) if a.parent.name in sym}

    needed = set()
    for e in ops:
        needed |= used(e)
    for name, _ in reversed(arrays):
        for i, e in enumerate(exprs[name]):
            if sym[name][i] in needed:
                needed |= used(e)
    keep = {name: [i for i in range(len(exprs[name])) if sym[name][i] in needed] for name, _ in arrays}
    remap, new = {}, {}
    for name, idxs in keep.items():
        if idxs:
            new[name] = sp.MatrixSymbol(name, len(idxs), 1)
            remap.update({sym[name][i]: new[name][k] for k, i in enumerate(idxs)})
    out = [(name, sp.Matrix([exprs[name][i].xreplace(remap) for i in keep[name]])) for name, _ in arrays
           if keep[name]]
    return out, [e.xreplace(remap) for e in ops]


def _wrappers_c(headers, orders):
    """Order-dispatch wrapper declarations and definitions (C/C++)."""
    from .writer import _parameter_names

    decls, defs = [], []
    first = f"_{min(orders)}"
    for decl in (h for h in headers if h.strip()):
        base = decl.split("(")[0].split()[-1]
        if not base.endswith(first) or base.endswith("batch" + first):
            continue
        stem = base[: -len(first)]
        wd = decl.replace(")", ", int order)").replace(base, stem)
        args = ", ".join(_parameter_names(decl))
        code = wd + " {\n  switch (order) {\n"
        for i in orders:
            code += f"  case {i}:\n    {stem}_{i}({args});\n    break;\n"
        code += "  }\n}\n"
        decls.append(wd + ";")
        defs.append(code)
    return decls, defs


def generate_kernel_code(
    kernel,
    order,
    name,
    orders=None,
    language="c",
    CSE=True,
    potential=True,
    field=True,
    minpow=None,
    include_dir=None,
    src_dir=None,
    operators=("S2M", "M2M", "M2L", "L2L", "L2P", "M2P", "P2P"),
    report=None,
    source_order=0,
    nchannels=1,
    prefix="",
    planar=False,
):
    """Write {name}.h/.c (or .cpp) -- or {name}.f90 -- for `kernel`.

    Emits, for every p in `orders` (default source_order+1..order): S2M_p,
    M2M_p, M2L_p, L2L_p, L2P_p, M2P_p; plus P2P and P2P_batch, and order-dispatch
    wrappers. Argument lists match the Laplace operators (dimension permitting):
        S2M_p(x,y,z, S, M)         M2M_p(x,y,z, M, Ms)    M2L_p(x,y,z, M, L)
        L2L_p(x,y,z, L, Ls)        L2P_p(x,y,z, L, F)     M2P_p(x,y,z, M, F)
        P2P(x,y,z, S, F)           P2P_batch(tx,ty,tz, sx,sy,sz, S, begin, end, F)
    F holds [phi, -dphi/dx_i...] (potential and/or field).

    `order` is the highest order emitted. FMMGEN_MAXORDER in the header is one more, an exclusive bound as in
    `generate_code`, so a driver loops `for (order = FMMGEN_MINORDER; order < FMMGEN_MAXORDER; ...)`.

    minpow: integer powers up to this are written as products instead of pow()
    calls, which is much faster for the steep kernels (r^-6, r^-12). None (default)
    uses 64, enough for every library kernel at the orders that are practical to
    generate; 0 keeps every pow() call.

    source_order s (0 charges, 1 dipoles, 2 quadrupoles, ...): a source has
    shell = C(s+d-1, d-1) coefficients S_m, |m| = s, and potential
    sum_m S_m d^m G(x - y). The multipole array holds s <= |n| <= p and the local
    array |n| <= p - s, exactly as in the Laplace source_order path.

    nchannels C > 1 stacks C independent source sets (S, M, L, F channel-major)
    that share one derivative tensor.

    prefix: prepended to every emitted function and to the header macros (FMMGEN_<PREFIX>MINORDER ...),
    so that operator sets for different kernels can be linked into one program. C and C++ only: Fortran
    output is namespaced by its module, which is named after `name`.

    planar: also emit the planar variants S2Mxy_p ... M2Pxy_p, P2Pxy and P2P_batchxy for sources and targets
    confined to the z = 0 plane (three dimensions; moments stay full 3D), valid for any kernel. They take the
    same arguments but smaller multipole and local arrays (n_z <= source_order, and n_z <= 1 with the
    field, n_z <= 0 without): see FMMGEN_PLANAR_MULTIPOLESIZE and FMMGEN_PLANAR_LOCALSIZE, indexed by order.

    report: optional dict, filled with {(op, p): (opscount, seconds)}.
    """
    if planar and kernel.dim != 3:
        raise ValueError("planar variants need a three-dimensional kernel")
    if prefix and language == "fortran":
        raise ValueError("prefix is not supported for fortran output (its module name is the namespace)")
    if prefix and not prefix.replace("_", "a").isalnum():
        raise ValueError("prefix must be a valid C identifier fragment")
    if language not in ("c", "c++", "fortran"):
        raise ValueError("language must be c, c++ or fortran")
    s, nch = source_order, nchannels
    dim, vars_ = kernel.dim, kernel.variables
    orders = list(orders) if orders is not None else list(range(s + 1, order + 1))
    if min(orders) < s:
        raise ValueError("every order must be >= source_order")
    order = max(orders)
    start = min(orders)
    p = FunctionPrinter(language=language, debug=not CSE, minpow=64 if minpow is None else minpow)
    coords = list(vars_)
    header, body = "", ""
    if report is None:
        report = {}
    shell = len(shell_indices(dim, s))

    def emit(opname, lhs, exprs, inputs, operator="+=", internal=None, **kw):
        nonlocal header, body
        t0 = time.time()
        head, code, ops = p.generate(prefix + opname, lhs, sp.Matrix(exprs), inputs, operator=operator,
                                     internal=internal or [], horner=False)
        if language != "fortran":
            code = _dedupe_R(code)
        header += head
        body += code + "\n"
        report[tuple(opname.rsplit("_", 1)) if opname[-1].isdigit() else (opname, None)] = (
            ops, time.time() - t0)

    zsym = vars_[-1]

    for q in orders:
        nm = len(source_indices(dim, q, s))
        nl = nterms(dim, q - s)
        Mq = sp.MatrixSymbol("M", nm * nch, 1)
        Lq = sp.MatrixSymbol("L", nl * nch, 1)
        Ssym = sp.MatrixSymbol("S", shell * nch, 1)
        if planar:
            keep_m, keep_l = _planar_keep(dim, q, s, nch, field)
            Mxy, Lxy = sp.MatrixSymbol("M", len(keep_m), 1), sp.MatrixSymbol("L", len(keep_l), 1)
            all_rows = None  # outputs F keep every row
            in_m, in_l = {"M": (nm * nch, keep_m)}, {"L": (nl * nch, keep_l)}

        def planar_emit(opname, lhs, ops, rows, inputs_map, inputs, Dt=None):
            ops_xy = _planar_ops(ops, range(len(ops)) if rows is None else rows, inputs_map, zsym)
            if Dt is None:
                emit(f"{opname}xy_{q}", lhs, ops_xy, inputs)
            else:
                arrays, ops_xy = _planar_arrays(Dt.arrays, ops_xy, zsym)
                emit(f"{opname}xy_{q}", lhs, ops_xy, inputs, internal=arrays)

        if "S2M" in operators:
            ops = _ops_S2M(dim, vars_, q, s, nch)
            emit(f"S2M_{q}", "M", ops, coords + [Ssym])
            if planar:
                planar_emit("S2M", "M", ops, keep_m, {}, coords + [Ssym])
        if "M2M" in operators:
            ops = _ops_M2M(dim, vars_, q, s, nch)
            emit(f"M2M_{q}", "Ms", ops, coords + [Mq])
            if planar:
                planar_emit("M2M", "Ms", ops, keep_m, in_m, coords + [Mxy])
        if "M2L" in operators:
            t0 = time.time()
            Dt = _DTensor(kernel, q)
            tD = time.time() - t0
            ops = _ops_M2L(dim, q, Dt, s, nch)
            emit(f"M2L_{q}", "L", ops, coords + [Mq], internal=Dt.arrays)
            report[("D-derivation", str(q))] = (0, tD)
            if planar:
                planar_emit("M2L", "L", ops, keep_l, in_m, coords + [Mxy], Dt)
        if "L2L" in operators:
            ops = _ops_L2L(dim, vars_, q, s, nch)
            emit(f"L2L_{q}", "Ls", ops, coords + [Lq])
            if planar:
                planar_emit("L2L", "Ls", ops, keep_l, in_l, coords + [Lxy])
        if "L2P" in operators:
            ops = _ops_L2P(dim, vars_, q, potential, field, s, nch)
            emit(f"L2P_{q}", "F", ops, coords + [Lq])
            if planar:
                planar_emit("L2P", "F", ops, all_rows, in_l, coords + [Lxy])
        if "M2P" in operators:
            Dt = _DTensor(kernel, q + (1 if field else 0))
            ops = _ops_M2P(dim, q, Dt, potential, field, s, nch)
            emit(f"M2P_{q}", "F", ops, coords + [Mq], internal=Dt.arrays)
            if planar:
                planar_emit("M2P", "F", ops, all_rows, in_m, coords + [Mxy], Dt)

    if "P2P" in operators:
        P2P = _ops_P2P(kernel, potential, field, s, nch)
        emit("P2P", "F", P2P, coords + [sp.MatrixSymbol("S", shell * nch, 1)])
        head, code, _ = p.generate_batch(prefix + "P2P_batch", "F", sp.Matrix(P2P), [str(v) for v in vars_],
                                         shell * nch)
        header += head
        body += _dedupe_R(code) + "\n" if language != "fortran" else code + "\n"
        if planar:
            P2P_xy = _planar_ops(P2P, range(len(P2P)), {}, zsym)
            emit("P2Pxy", "F", P2P_xy, coords + [sp.MatrixSymbol("S", shell * nch, 1)])
            head, code, _ = p.generate_batch(prefix + "P2P_batchxy", "F", sp.Matrix(P2P_xy),
                                             [str(v) for v in vars_], shell * nch)
            header += head
            body += _dedupe_R(code) + "\n" if language != "fortran" else code + "\n"

    osize = ((1 if potential else 0) + (dim if field else 0)) * nch
    planar_tables = None
    if planar:  # indexed by order 0..order; 0 where the order is below the source order
        sizes = [planar_sizes(dim, o, s, nch, field) if o >= s else (0, 0) for o in range(order + 1)]
        planar_tables = ([a for a, _ in sizes], [b for _, b in sizes])

    if language == "fortran":
        _write_fortran_kernel(p, body, name, orders, kernel, osize, src_dir or include_dir,
                              planar_sizes=planar_tables if planar else None)
        return report

    decls, defs = _wrappers_c(header.split(";\n"), orders)
    header += "\n".join(decls) + "\n"
    body += "\n".join(defs)
    hext = "h"
    fext = "c" if language == "c" else "cpp"
    from .writer import _RESTRICT_MACRO

    def path(d, ext):
        return f"{d.rstrip('/')}/{name}.{ext}" if d else f"{name}.{ext}"

    with open(path(include_dir, hext), "w") as f:
        f.write("#pragma once\n")
        f.write("#include <cstddef>\n" if language == "c++" else "#include <stddef.h>\n")
        f.write(_RESTRICT_MACRO)
        m = "FMMGEN_" + (prefix.upper().rstrip("_") + "_" if prefix else "")
        f.write(f"#define {m}MINORDER {start}\n#define {m}MAXORDER {order + 1}\n")
        f.write(f"#define {m}DIM {dim}\n#define {m}SOURCEORDER {s}\n#define {m}NCHANNELS {nch}\n")
        f.write(f"#define {m}SOURCESIZE {shell * nch}\n#define {m}OUTPUTSIZE {osize}\n")
        if planar:
            f.write(f"#define {m}PLANAR 1\n")
            f.write("/* Array sizes of the xy-suffixed operators only, indexed by order. */\n")
            f.write(f"static const size_t {m}PLANAR_MULTIPOLESIZE[] = {{{', '.join(map(str, planar_tables[0]))}}};\n")
            f.write(f"static const size_t {m}PLANAR_LOCALSIZE[] = {{{', '.join(map(str, planar_tables[1]))}}};\n")
        f.write(f"/* kernel: {kernel.name.replace('*/', '')} ; derivative method: {kernel.method} */\n")
        f.write(header)
    with open(path(src_dir, fext), "w") as f:
        f.write(f'#include "{name}.{hext}"\n')
        f.write('#include "math.h"\n' if language == "c" else "#include <cmath>\n")
        f.write(body)
    return report


def _write_fortran_kernel(p, body, name, orders, kernel, osize, out_dir, planar_sizes=None):
    order, start = max(orders), min(orders)
    out = [f"module {name}", "  implicit none", "  private", "  integer, parameter, public :: wp = kind(1.0d0)"]
    out += [
        f"  integer, parameter, public :: FMMGEN_MINORDER = {start}",
        f"  integer, parameter, public :: FMMGEN_MAXORDER = {order + 1}",
        f"  integer, parameter, public :: FMMGEN_DIM = {kernel.dim}",
        "  integer, parameter, public :: FMMGEN_SOURCESIZE = 1",
        f"  integer, parameter, public :: FMMGEN_OUTPUTSIZE = {osize}",
        f"  logical, parameter, public :: FMMGEN_PLANAR = {'.true.' if planar_sizes else '.false.'}",
    ]
    if planar_sizes:
        for label, values in zip(("MULTIPOLE", "LOCAL"), planar_sizes):
            out.append(f"  integer, parameter, public :: FMMGEN_PLANAR_{label}SIZE(0:{order}) = (/ "
                       + ", ".join(map(str, values)) + " /)")
    suffix = f"_{start}"
    stems = [n[: -len(suffix)] for n in p.signatures if n.endswith(suffix)]
    public = list(dict.fromkeys(stems + [n for n in p.signatures if not n.endswith(suffix) and not n[-1].isdigit()]))
    for k in range(0, len(public), 8):
        out.append("  public :: " + ", ".join(public[k:k + 8]))
    out += ["contains", "", body]
    for stem in stems:
        args = p.signatures[stem + suffix] + [("order", "int")]
        out.append(f"subroutine {stem}({', '.join(n for n, _ in args)})")
        out += ["  " + d for d in p._fortran_declarations(args)]
        out.append("  select case (order)")
        for i in orders:
            out.append(f"  case ({i})")
            out.append(f"    call {stem}_{i}({', '.join(n for n, _ in args[:-1])})")
        out += ["  case default", f"    stop 'fmmgen: {stem} called with an order that was not generated'",
                "  end select", f"end subroutine {stem}", ""]
    out.append(f"end module {name}")
    path = f"{out_dir.rstrip('/')}/{name}.f90" if out_dir else f"{name}.f90"
    with open(path, "w") as f:
        f.write("\n".join(out) + "\n")


# ----------------------------------------------------------------------------
# Kernel library: radial kernels common in simulation (all via Kernel.radial,
# so they all use the radial recurrence; none besides Laplace-type has a pde)
# ----------------------------------------------------------------------------
_r = sp.Symbol("r", positive=True)


def power_law(a, dim=3, coeff=1, **kw):
    """coeff * r^(-a): a=1 Coulomb/gravity, a=3 dipolar, a=6 dispersion, a=12 repulsion."""
    pde = None
    if a == dim - 2 and dim >= 3:
        pde = {tuple(2 if j == i else 0 for j in range(dim)): 1 for i in range(dim)}
    return Kernel.radial(coeff * _r ** (-sp.sympify(a)), _r, dim, pde=kw.pop("pde", pde), **kw)


def lennard_jones(epsilon=1.0, sigma=1.0, dim=3, **kw):
    """4 eps ((sigma/r)^12 - (sigma/r)^6): a sum of two power laws."""
    e, sg = _as_number(epsilon), _as_number(sigma)
    return Kernel.radial(4 * e * ((sg / _r) ** 12 - (sg / _r) ** 6), _r, dim, **kw)


def morse(D=1.0, alpha=1.0, r0=1.0, dim=3, **kw):
    """D (exp(-2 alpha (r - r0)) - 2 exp(-alpha (r - r0)))."""
    D, al, r0 = _as_number(D), _as_number(alpha), _as_number(r0)
    return Kernel.radial(D * (sp.exp(-2 * al * (_r - r0)) - 2 * sp.exp(-al * (_r - r0))), _r, dim, **kw)


def buckingham(A=1.0, rho=1.0, C=1.0, dim=3, **kw):
    """A exp(-r / rho) - C / r^6."""
    A, rho, C = _as_number(A), _as_number(rho), _as_number(C)
    return Kernel.radial(A * sp.exp(-_r / rho) - C / _r**6, _r, dim, **kw)


def ewald_real(alpha=1.0, dim=3, **kw):
    """erfc(alpha r) / r: the real-space part of an Ewald split of 1/r."""
    return Kernel.radial(sp.erfc(_as_number(alpha) * _r) / _r, _r, dim, **kw)


def gaussian_screened_coulomb(alpha=1.0, dim=3, **kw):
    """erf(alpha r) / r: Coulomb between a point charge and a Gaussian charge."""
    a = _as_number(alpha)
    return Kernel.radial(sp.erf(a * _r) / _r, _r, dim,
                         radial_derivs=lambda kmax: _stable_erf_derivs(a, kmax), **kw)


def softened_power(eps=0.1, a=1, dim=3, **kw):
    """(r^2 + eps^2)^(-a/2): Plummer-softened gravity for a=1."""
    return Kernel.radial((_r**2 + _as_number(eps) ** 2) ** (-sp.Rational(a, 2)), _r, dim, **kw)


def yukawa(k=1.0, dim=3, **kw):
    """exp(-k r) / r, with its Helmholtz-type pde."""
    k = _as_number(k)
    pde = {**{tuple(2 if j == i else 0 for j in range(dim)): 1 for i in range(dim)}, tuple([0] * dim): -k**2}
    return Kernel.radial(sp.exp(-k * _r) / _r, _r, dim, pde=kw.pop("pde", pde), **kw)


def _helmholtz_pde(k, dim):
    return {**{tuple(2 if j == i else 0 for j in range(dim)): 1 for i in range(dim)}, tuple([0] * dim): k**2}


def helmholtz_cos(k=1.0, dim=3, **kw):
    """cos(k r) / r: the real part of exp(i k r) / r, with its Helmholtz pde.

    Low-frequency Helmholtz only (k D up to about 4 across a box): the Cartesian
    expansion loses accuracy as the box grows in wavelengths. For a complex source
    q = a + i b, Re/Im of the potential are
        Re = cos * a - sin * b,   Im = cos * b + sin * a
    so generate each of `helmholtz_cos` and `helmholtz_sin` with nchannels=2
    (channels a and b), which shares each derivative tensor between the channels.
    """
    k = _as_number(k)
    return Kernel.radial(sp.cos(k * _r) / _r, _r, dim, pde=kw.pop("pde", _helmholtz_pde(k, dim)), **kw)


def helmholtz_sin(k=1.0, dim=3, **kw):
    """sin(k r) / r: the imaginary part of exp(i k r) / r (smooth at r = 0)."""
    k = _as_number(k)
    return Kernel.radial(sp.sin(k * _r) / _r, _r, dim, pde=kw.pop("pde", _helmholtz_pde(k, dim)), **kw)
