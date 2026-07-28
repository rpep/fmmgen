"""
Stage 2: compressed operators against the uncompressed path.

Each compressed operator must be numerically equivalent to the uncompressed
one under the appropriate translation of its input and output arrays:

    multipole array   compressed <-> full   by Dec^T (out) / lift-with-zeros (in)
    local array       compressed <-> full   by restriction (out) / Dec (in)

Equivalence is checked on random inputs rather than symbolically, because the
compressed forms are deliberately *not* term-by-term rearrangements of the
uncompressed ones -- they contract over a smaller index set.

Pure Python; no compiler needed.
"""
import numpy as np
import pytest
import sympy as sp

from fmmgen import generator as g
from fmmgen.harmonic import Nkeep, decompress, keep_mappings
from fmmgen.utils import Nterms, generate_mappings

symbols = sp.symbols("x y z")
x, y, z = symbols

ORDERS = [1, 2, 3, 4, 5]
TOL = 1e-9


def _setup(p):
    full, _ = generate_mappings(p, symbols, "grevlex", source_order=0)
    keep, _ = keep_mappings(p, symbols)
    return full, keep, decompress(p, symbols)


def _coords(rng):
    return {s: float(v) for s, v in zip(symbols, rng.standard_normal(3) * 0.3)}


def _trace_free_D(p, rng, full, keep, dec):
    """Random D satisfying D_{k+(0,0,2)} = -D_{k+(2,0,0)} - D_{k+(0,2,0)}."""
    free = {n: rng.standard_normal() for n in keep}
    vals = {n: sum(c * free[k] for k, c in dec[n].items()) for n in full}
    D = sp.MatrixSymbol("D", Nterms(p), 1)
    return {D[full[n]]: vals[n] for n in full}


def _evaluate(exprs, subs):
    return np.array([float(sp.N(e.xreplace(subs))) for e in exprs])


def _project(vals, full, keep, dec):
    """Dec^T applied to a full-indexed numeric vector."""
    out = {i: 0.0 for i in keep}
    for n, idx in full.items():
        for k, c in dec[n].items():
            out[k] += float(c) * vals[idx]
    return np.array([out[i] for i in keep])


def _restrict(vals, full, keep):
    return np.array([vals[full[i]] for i in keep])


def _lift_subs(name, p, full, keep, comp_vals):
    """Full-indexed substitution putting compressed values in retained slots."""
    F = sp.MatrixSymbol(name, Nterms(p), 1)
    return {
        F[idx]: (float(comp_vals[keep[n]]) if n in keep else 0.0)
        for n, idx in full.items()
    }


def _expand_subs(name, p, full, keep, dec, comp_vals):
    """Full-indexed substitution reconstructing a local array from Dec."""
    F = sp.MatrixSymbol(name, Nterms(p), 1)
    return {
        F[idx]: float(sum(c * comp_vals[keep[k]] for k, c in dec[n].items()))
        for n, idx in full.items()
    }


def _comp_subs(name, p, comp_vals):
    C = sp.MatrixSymbol(name, Nkeep(p), 1)
    return {C[i]: float(v) for i, v in enumerate(comp_vals)}


# --------------------------------------------------------------------------

@pytest.mark.parametrize("p", ORDERS)
def test_s2m(p):
    """Compressed S2M == Dec^T of uncompressed S2M."""
    rng = np.random.default_rng(100 + p)
    full, keep, dec = _setup(p)
    coords = _coords(rng)

    S = sp.MatrixSymbol("S", 1, 1)
    s_subs = {S[0]: float(rng.standard_normal())}

    ref = _evaluate(g.generate_S2M_operators(p, symbols, full), {**coords, **s_subs})
    got = _evaluate(
        g.generate_S2M_operators_compressed(p, symbols, full, keep),
        {**coords, **s_subs},
    )
    assert np.allclose(got, _project(ref, full, keep, dec), atol=TOL)


@pytest.mark.parametrize("p", ORDERS)
def test_m2m(p):
    """Compressed M2M agrees after lifting the input and projecting the output."""
    rng = np.random.default_rng(200 + p)
    full, keep, dec = _setup(p)
    coords = _coords(rng)
    Mt = rng.standard_normal(Nkeep(p))

    ref = _evaluate(
        g.generate_M_shift_operators(p, symbols, full),
        {**coords, **_lift_subs("M", p, full, keep, Mt)},
    )
    got = _evaluate(
        g.generate_M_shift_operators_compressed(p, symbols, full, keep),
        {**coords, **_comp_subs("M", p, Mt)},
    )
    assert np.allclose(got, _project(ref, full, keep, dec), atol=TOL)


@pytest.mark.parametrize("p", ORDERS)
def test_m2l(p):
    """Compressed M2L agrees after lifting the multipole and restricting the local."""
    rng = np.random.default_rng(300 + p)
    full, keep, dec = _setup(p)
    Mt = rng.standard_normal(Nkeep(p))
    d_subs = _trace_free_D(p, rng, full, keep, dec)

    ref = _evaluate(
        g.generate_L_operators(p, symbols, full, full),
        {**d_subs, **_lift_subs("M", p, full, keep, Mt)},
    )
    got = _evaluate(
        g.generate_L_operators_compressed(p, symbols, full, full, keep),
        {**d_subs, **_comp_subs("M", p, Mt)},
    )
    assert np.allclose(got, _restrict(ref, full, keep), atol=TOL)

    # and the restricted local must decompress back to the full one, which is
    # the property that lets L2L and L2P consume the compressed array
    rebuilt = np.array(
        [float(sum(c * got[keep[k]] for k, c in dec[n].items())) for n in full]
    )
    assert np.allclose(rebuilt, ref, atol=TOL)


@pytest.mark.parametrize("p", ORDERS)
def test_l2l(p):
    """Compressed L2L agrees after expanding the input and restricting the output."""
    rng = np.random.default_rng(400 + p)
    full, keep, dec = _setup(p)
    coords = _coords(rng)
    Lt = rng.standard_normal(Nkeep(p))

    ref = _evaluate(
        g.generate_L_shift_operators(p, symbols, full),
        {**coords, **_expand_subs("L", p, full, keep, dec, Lt)},
    )
    got = _evaluate(
        g.generate_L_shift_operators_compressed(p, symbols, full, keep),
        {**coords, **_comp_subs("L", p, Lt)},
    )
    assert np.allclose(got, _restrict(ref, full, keep), atol=TOL)


@pytest.mark.parametrize("p", ORDERS)
def test_l2p(p):
    """Potential and field from a compressed local expansion are unchanged."""
    rng = np.random.default_rng(500 + p)
    full, keep, dec = _setup(p)
    coords = _coords(rng)
    Lt = rng.standard_normal(Nkeep(p))

    ref = _evaluate(
        g.generate_L2P_operators(p, symbols, full),
        {**coords, **_expand_subs("L", p, full, keep, dec, Lt)},
    )
    got = _evaluate(
        g.generate_L2P_operators_compressed(p, symbols, full, keep),
        {**coords, **_comp_subs("L", p, Lt)},
    )
    assert np.allclose(got, ref, atol=TOL)


@pytest.mark.parametrize("p", ORDERS)
def test_m2p(p):
    """Potential and field from a compressed multipole are unchanged."""
    rng = np.random.default_rng(600 + p)
    full, keep, dec = _setup(p)
    # M2P evaluates 1/R derivatives inline, so keep the target well separated.
    # Phi_derivatives emits the atomic symbol Rinv rather than 1/sqrt(...), so
    # it has to be supplied consistently with the coordinates.
    coords = {x: 1.7, y: -2.3, z: 0.9}
    coords[sp.Symbol("Rinv")] = 1.0 / np.sqrt(1.7**2 + 2.3**2 + 0.9**2)
    Mt = rng.standard_normal(Nkeep(p))

    ref = _evaluate(
        g.generate_M2P_operators(p, symbols, full),
        {**coords, **_lift_subs("M", p, full, keep, Mt)},
    )
    got = _evaluate(
        g.generate_M2P_operators_compressed(p, symbols, full, keep),
        {**coords, **_comp_subs("M", p, Mt)},
    )
    assert np.allclose(got, ref, atol=TOL)


def test_source_order_is_rejected():
    """The compressed path is source_order = 0 only, and says so."""
    full, keep, _ = _setup(3)
    with pytest.raises(NotImplementedError, match="source_order"):
        g.generate_M_shift_operators_compressed(3, symbols, full, keep, source_order=1)
