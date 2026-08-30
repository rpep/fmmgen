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
from fmmgen.harmonic import keep_mappings
from fmmgen.utils import Nterms, generate_mappings

symbols = sp.symbols("x y z")
x, y, z = symbols

ORDERS = [1, 2, 3, 4, 5]
TOL = 1e-9

# source_order 0 (charges) and 1 (point dipoles, the micromagnetics case).
# At s > 0 the multipole spans s <= |n| <= p while the local spans |m| <= p-s,
# so the two arrays have different retained sets -- which is the whole reason
# these are threaded separately.
CASES = [(p, s) for s in (0, 1) for p in ORDERS if p > s]


def _setup(p, s=0):
    M_dict, _ = generate_mappings(p, symbols, "grevlex", source_order=s)
    L_dict, _ = generate_mappings(p - s, symbols, "grevlex", source_order=0)
    keep_M, _ = keep_mappings(p, symbols, source_order=s)
    keep_L, _ = keep_mappings(p - s, symbols, source_order=0)
    from fmmgen.harmonic import decompress_for
    return M_dict, L_dict, keep_M, keep_L, decompress_for(M_dict), decompress_for(L_dict)


def _coords(rng):
    return {s: float(v) for s, v in zip(symbols, rng.standard_normal(3) * 0.3)}


def _trace_free_D(p, rng, M_dict, keep_M, dec_M):
    """
    Random D satisfying D_{k+(0,0,2)} = -D_{k+(2,0,0)} - D_{k+(0,2,0)}.

    Indexed by M_dict, not by the full 0 <= |n| <= p set: expansions.L reads
    the derivative array as D[M_dict[n+m]], so at source_order > 0 the slots
    follow the multipole index map. Legitimate because the relation preserves
    total degree, so an entry of degree >= s only ever refers to others of the
    same degree, which are also in M_dict.
    """
    free = {n: rng.standard_normal() for n in keep_M}
    vals = {n: sum(c * free[k] for k, c in dec_M[n].items()) for n in M_dict}
    D = sp.MatrixSymbol("D", Nterms(p), 1)
    return {D[M_dict[n]]: vals[n] for n in M_dict}


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


def _comp_subs(name, comp_vals):
    C = sp.MatrixSymbol(name, len(comp_vals), 1)
    return {C[i]: float(v) for i, v in enumerate(comp_vals)}


# --------------------------------------------------------------------------

@pytest.mark.parametrize("p,s", CASES)
def test_s2m(p, s):
    """Compressed S2M == Dec^T of uncompressed S2M."""
    rng = np.random.default_rng(100 + p + 50 * s)
    full, L_dict, keep_M, keep_L, dec_M, dec_L = _setup(p, s)
    coords = _coords(rng)

    nsrc = len([n for n in full if sum(n) == s])
    S = sp.MatrixSymbol("S", nsrc, 1)
    s_subs = {S[i]: float(v) for i, v in enumerate(rng.standard_normal(nsrc))}

    ref = _evaluate(g.generate_S2M_operators(p, symbols, full, source_order=s),
                    {**coords, **s_subs})
    got = _evaluate(
        g.generate_S2M_operators_compressed(p, symbols, full, keep_M, source_order=s),
        {**coords, **s_subs},
    )
    assert np.allclose(got, _project(ref, full, keep_M, dec_M), atol=TOL)


@pytest.mark.parametrize("p,s", CASES)
def test_m2m(p, s):
    """Compressed M2M agrees after lifting the input and projecting the output."""
    rng = np.random.default_rng(200 + p + 50 * s)
    full, L_dict, keep_M, keep_L, dec_M, dec_L = _setup(p, s)
    coords = _coords(rng)
    Mt = rng.standard_normal(len(keep_M))

    ref = _evaluate(
        g.generate_M_shift_operators(p, symbols, full, source_order=s),
        {**coords, **_lift_subs("M", p, full, keep_M, Mt)},
    )
    got = _evaluate(
        g.generate_M_shift_operators_compressed(p, symbols, full, keep_M, source_order=s),
        {**coords, **_comp_subs("M", Mt)},
    )
    assert np.allclose(got, _project(ref, full, keep_M, dec_M), atol=TOL)


@pytest.mark.parametrize("p,s", CASES)
def test_m2l(p, s):
    """Compressed M2L agrees after lifting the multipole and restricting the local."""
    rng = np.random.default_rng(300 + p + 50 * s)
    full, L_dict, keep_M, keep_L, dec_M, dec_L = _setup(p, s)
    Mt = rng.standard_normal(len(keep_M))
    d_subs = _trace_free_D(p, rng, full, keep_M, dec_M)

    ref = _evaluate(
        g.generate_L_operators(p, symbols, full, L_dict, source_order=s),
        {**d_subs, **_lift_subs("M", p, full, keep_M, Mt)},
    )
    got = _evaluate(
        g.generate_L_operators_compressed(p, symbols, full, L_dict, keep_M, keep_L,
                                          source_order=s),
        {**d_subs, **_comp_subs("M", Mt)},
    )
    assert np.allclose(got, _restrict(ref, L_dict, keep_L), atol=TOL)

    # and the restricted local must decompress back to the full one, which is
    # the property that lets L2L and L2P consume the compressed array
    rebuilt = np.array(
        [float(sum(c * got[keep_L[k]] for k, c in dec_L[n].items())) for n in L_dict]
    )
    assert np.allclose(rebuilt, ref, atol=TOL)


@pytest.mark.parametrize("p,s", CASES)
def test_l2l(p, s):
    """Compressed L2L agrees after expanding the input and restricting the output."""
    rng = np.random.default_rng(400 + p + 50 * s)
    full, L_dict, keep_M, keep_L, dec_M, dec_L = _setup(p, s)
    coords = _coords(rng)
    Lt = rng.standard_normal(len(keep_L))

    ref = _evaluate(
        g.generate_L_shift_operators(p, symbols, L_dict, source_order=s),
        {**coords, **_expand_subs("L", p, L_dict, keep_L, dec_L, Lt)},
    )
    got = _evaluate(
        g.generate_L_shift_operators_compressed(p, symbols, L_dict, keep_L, source_order=s),
        {**coords, **_comp_subs("L", Lt)},
    )
    assert np.allclose(got, _restrict(ref, L_dict, keep_L), atol=TOL)


@pytest.mark.parametrize("p,s", CASES)
def test_l2p(p, s):
    """Potential and field from a compressed local expansion are unchanged."""
    rng = np.random.default_rng(500 + p + 50 * s)
    full, L_dict, keep_M, keep_L, dec_M, dec_L = _setup(p, s)
    coords = _coords(rng)
    Lt = rng.standard_normal(len(keep_L))

    ref = _evaluate(
        g.generate_L2P_operators(p, symbols, L_dict),
        {**coords, **_expand_subs("L", p, L_dict, keep_L, dec_L, Lt)},
    )
    got = _evaluate(
        g.generate_L2P_operators_compressed(p, symbols, L_dict, keep_L),
        {**coords, **_comp_subs("L", Lt)},
    )
    assert np.allclose(got, ref, atol=TOL)


@pytest.mark.parametrize("p,s", CASES)
def test_m2p(p, s):
    """Potential and field from a compressed multipole are unchanged."""
    rng = np.random.default_rng(600 + p + 50 * s)
    full, L_dict, keep_M, keep_L, dec_M, dec_L = _setup(p, s)
    # M2P evaluates 1/R derivatives inline, so keep the target well separated.
    # Phi_derivatives emits the atomic symbol Rinv rather than 1/sqrt(...), so
    # it has to be supplied consistently with the coordinates.
    coords = {x: 1.7, y: -2.3, z: 0.9}
    coords[sp.Symbol("Rinv")] = 1.0 / np.sqrt(1.7**2 + 2.3**2 + 0.9**2)
    Mt = rng.standard_normal(len(keep_M))

    ref = _evaluate(
        g.generate_M2P_operators(p, symbols, full, source_order=s),
        {**coords, **_lift_subs("M", p, full, keep_M, Mt)},
    )
    got = _evaluate(
        g.generate_M2P_operators_compressed(p, symbols, full, keep_M, source_order=s),
        {**coords, **_comp_subs("M", Mt)},
    )
    assert np.allclose(got, ref, atol=TOL)


def test_quadrupole_sources_are_rejected():
    """
    source_order >= 2 is refused rather than silently miscomputed: at degree 2
    the source shell itself compresses (6 monomials -> 5 retained), so S2M
    would have to project the caller's moment array.
    """
    full, _, keep_M, _, _, _ = _setup(4, 2)
    with pytest.raises(NotImplementedError, match="source_order"):
        g.generate_M_shift_operators_compressed(4, symbols, full, keep_M, source_order=2)
