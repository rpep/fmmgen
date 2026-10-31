"""The derivative-tensor evaluators (recurrence, Faa di Bruno, symbolic, and
the pde= recurrences on top) against sp.diff -- pure Python, no compiler."""
import math

import pytest
import sympy as sp
from sympy.matrices.expressions.matexpr import MatrixElement

from fmmgen.kernel import Kernel, _DTensor, multi_indices
from kernel_util import KERNELS, LAP3, r3


def evaluate_D(kernel, p, pt):
    """Evaluate every D_n by executing the generated array definitions."""
    Dt = _DTensor(kernel, p)
    R = math.sqrt(sum(v * v for v in pt))
    vals = {}
    sub = {v: sp.Float(c) for v, c in zip(kernel.variables, pt)}
    sub[sp.Symbol("R")] = sp.Float(R)
    for name, mat in Dt.arrays:
        out = []
        for e in mat:
            e = sp.sympify(e).xreplace(sub)
            e = sp.sympify(e)
            e = e.xreplace({a: sp.Float(vals[(a.parent.name, int(a.i))]) for a in e.atoms(MatrixElement)})
            out.append(float(e))
        for i, v in enumerate(out):
            vals[(name, i)] = v
    return {n: vals[(loc.parent.name, int(loc.i))] for n, loc in Dt.loc.items()}


def reference(kernel, p, pt):
    sub = dict(zip(kernel.variables, pt))
    out = {}
    for n in multi_indices(kernel.dim, p):
        d = kernel.expr
        for v, e in zip(kernel.variables, n):
            if e:
                d = sp.diff(d, v, e)
        out[n] = float(d.subs(sub))
    return out


CASES = [
    ("laplace", 3, [0.7, -0.4, 0.9]),
    ("yukawa", 3, [0.7, -0.4, 0.9]),
    ("r", 3, [0.7, -0.4, 0.9]),
    ("dipolelike", 3, [0.7, -0.4, 0.9]),
    ("plummer", 3, [0.7, -0.4, 0.9]),
    ("gaussian", 3, [0.7, -0.4, 0.9]),
    ("log2d", 2, [0.7, -0.4]),
]


@pytest.mark.parametrize("method", ["recurrence", "faa", "symbolic"])
@pytest.mark.parametrize("name,dim,pt", CASES)
def test_D_matches_sympy(name, dim, pt, method):
    k = KERNELS[name](method=method)
    p = 4
    got, ref = evaluate_D(k, p, pt), reference(k, p, pt)
    for n in ref:
        assert got[n] == pytest.approx(ref[n], rel=1e-9, abs=1e-10), (name, method, n)


@pytest.mark.parametrize("name,dim,pt", CASES[:3] + CASES[-1:])
def test_D_with_pde_matches_sympy(name, dim, pt):
    """The pde= path must agree with the pde-free one (both vs sp.diff)."""
    k = KERNELS[name]()
    assert k.pde
    got, ref = evaluate_D(k, 5, pt), reference(k, 5, pt)
    for n in ref:
        assert got[n] == pytest.approx(ref[n], rel=1e-8, abs=1e-9), (name, n)


def test_wrong_pde_rejected():
    with pytest.raises(ValueError):
        Kernel(sp.exp(-r3) / r3, pde=LAP3)  # not harmonic


def test_derivs_hook():
    """A user-supplied scheme: here D_n by sp.diff, in a single layer."""
    def scheme(kernel, pmax, idx):
        return [[(n, sp.diff(kernel.expr, *[a for v, e in zip(kernel.variables, n) for a in ((v, e) if e else ())])
                  if sum(n) else kernel.expr) for n in idx]]
    k = Kernel(1 / r3, derivs=scheme)
    got, ref = evaluate_D(k, 3, [0.7, -0.4, 0.9]), reference(k, 3, [0.7, -0.4, 0.9])
    for n in ref:
        assert got[n] == pytest.approx(ref[n], rel=1e-9)


def test_recurrence_table_size():
    """Table entries are C(p+d+1, d+1), not the exploding symbolic size."""
    k = Kernel(1 / r3)
    Dt = _DTensor(k, 8)
    entries = sum(len(m) for _, m in Dt.arrays)
    assert entries == math.comb(8 + 4, 4)
