"""
Stage 1 gate for harmonic (traceless) compression of Cartesian expansions.

The claim under test is that the M2L matrix H_{nm} = D_{n+m}, which is Hankel,
factors exactly through the retained index set:

    H = Dec * H[keep, keep] * Dec^T

If that holds, compressed M2L is literally the submatrix H[keep, keep] applied
to Mtilde = Dec^T M, and no separate construction is needed. If it fails --
most plausibly at the boundary of the |n| + |m| <= p truncation -- the whole
approach collapses back to the general construction of Fernando and Klockner
(arXiv:2305.17867) section 4.1.

These are pure Python and need no compiler.
"""
from math import comb

import numpy as np
import pytest
import sympy as sp

from fmmgen.harmonic import Nkeep, decompress, expand_index, keep_mappings
from fmmgen.utils import Nterms, generate_mappings

symbols = sp.symbols("x y z")

ORDERS = list(range(1, 13))


# --------------------------------------------------------------------------
# counting
# --------------------------------------------------------------------------

@pytest.mark.parametrize("p", ORDERS)
def test_retained_count_is_harmonic_dimension(p):
    """(p+1)^2 retained, matching C(p+3,3) - C(p+1,3) and the (2l+1) sum."""
    assert Nkeep(p) == (p + 1) ** 2
    assert Nkeep(p) == comb(p + 3, 3) - comb(p + 1, 3)
    assert Nkeep(p) == sum(2 * l + 1 for l in range(p + 1))

    fwd, rev = keep_mappings(p, symbols)
    assert len(fwd) == Nkeep(p)
    assert len(rev) == Nkeep(p)
    assert Nterms(p) == comb(p + 3, 3)


@pytest.mark.parametrize("p", ORDERS)
def test_keep_mappings_are_consistent(p):
    """Indices are contiguous from zero, retained-only, and order-preserving."""
    fwd, rev = keep_mappings(p, symbols)
    assert sorted(fwd.values()) == list(range(Nkeep(p)))
    assert all(rev[i] == n for n, i in fwd.items())
    assert all(n[2] <= 1 for n in fwd)
    assert all(sum(n) <= p for n in fwd)

    # relative ordering must match the uncompressed array, so that degree
    # shells stay contiguous and prefix-closed as elsewhere in fmmgen
    full, _ = generate_mappings(p, symbols)
    kept_in_full_order = [n for n in full if n[2] <= 1]
    assert list(fwd.keys()) == kept_in_full_order


# --------------------------------------------------------------------------
# the elimination itself
# --------------------------------------------------------------------------

def test_expand_index_known_cases():
    assert expand_index((1, 0, 1)) == {(1, 0, 1): 1}
    assert expand_index((0, 0, 2)) == {(2, 0, 0): -1, (0, 2, 0): -1}
    assert expand_index((0, 0, 3)) == {(2, 0, 1): -1, (0, 2, 1): -1}
    # d_z^4 = (d_x^2 + d_y^2)^2
    assert expand_index((0, 0, 4)) == {(4, 0, 0): 1, (2, 2, 0): 2, (0, 4, 0): 1}


@pytest.mark.parametrize("p", ORDERS)
def test_expansion_preserves_total_degree(p):
    """
    Every elimination stays within its degree shell.

    This is what makes compression compatible with the |n| + |m| <= p
    truncation: no index of higher degree is ever referenced, so the
    truncated support is closed under the relation.
    """
    for n, terms in decompress(p, symbols).items():
        assert all(sum(k) == sum(n) for k in terms), n
        assert all(k[2] <= 1 for k in terms), n


@pytest.mark.parametrize("p", [1, 2, 3, 4, 5, 6])
def test_expansion_matches_real_derivatives_of_one_over_r(p):
    """
    Check the relation against actual derivatives of 1/R, to pin the sign
    convention against sympy rather than against our own algebra.
    """
    x, y, z = symbols
    R = sp.sqrt(x**2 + y**2 + z**2)
    phi = 1 / R
    subs = {x: sp.Rational(3, 7), y: sp.Rational(-5, 11), z: sp.Rational(2, 3)}

    cache = {}

    def D(n):
        if n not in cache:
            cache[n] = sp.diff(phi, x, n[0], y, n[1], z, n[2]).subs(subs)
        return cache[n]

    for n, terms in decompress(p, symbols).items():
        lhs = D(n)
        rhs = sum(coeff * D(k) for k, coeff in terms.items())
        assert sp.simplify(lhs - rhs) == 0, n


# --------------------------------------------------------------------------
# the gate
# --------------------------------------------------------------------------

def _random_harmonic_derivatives(p, rng):
    """
    A random vector D indexed by |n| <= p satisfying the trace-free relation.

    Free values are assigned on the retained set and the rest filled in by the
    relation. The identity under test is linear in D, so a random valid D
    tests it generically -- far cheaper than differentiating 1/R to order 12.
    """
    keep_fwd, _ = keep_mappings(p, symbols)
    free = {n: rng.standard_normal() for n in keep_fwd}
    return {
        n: sum(c * free[k] for k, c in terms.items())
        for n, terms in decompress(p, symbols).items()
    }


@pytest.mark.parametrize("p", ORDERS)
def test_random_derivative_vector_is_trace_free(p):
    """Sanity check on the fixture before it is used to test anything else."""
    rng = np.random.default_rng(0)
    D = _random_harmonic_derivatives(p, rng)
    for n in decompress(p, symbols):
        if sum(n) + 2 > p:
            continue
        lap = (
            D[(n[0] + 2, n[1], n[2])]
            + D[(n[0], n[1] + 2, n[2])]
            + D[(n[0], n[1], n[2] + 2)]
        )
        assert abs(lap) < 1e-10, (n, lap)


@pytest.mark.parametrize("p", ORDERS)
def test_m2l_factors_through_the_retained_submatrix(p):
    """
    THE GATE.

        L_n = sum_{|m| <= p - |n|} M_m D_{n+m}

    computed directly, against the compressed route

        Mtilde = Dec^T M,  Ltilde = H[keep,keep] Mtilde,  L = Dec Ltilde.

    Failure here means H does not factor through the retained set under the
    truncation, and the submatrix shortcut is invalid.
    """
    rng = np.random.default_rng(p)
    full, _ = generate_mappings(p, symbols)
    keep_fwd, _ = keep_mappings(p, symbols)
    dec = decompress(p, symbols)
    D = _random_harmonic_derivatives(p, rng)

    M = {m: rng.standard_normal() for m in full}

    # direct, uncompressed
    L = {}
    for n in full:
        L[n] = sum(
            M[m] * D[(n[0] + m[0], n[1] + m[1], n[2] + m[2])]
            for m in full
            if sum(m) + sum(n) <= p
        )

    # compressed: project, contract on the retained submatrix, expand
    Mtilde = {j: 0.0 for j in keep_fwd}
    for m in full:
        for j, coeff in dec[m].items():
            Mtilde[j] += coeff * M[m]

    Ltilde = {}
    for i in keep_fwd:
        Ltilde[i] = sum(
            D[(i[0] + j[0], i[1] + j[1], i[2] + j[2])] * Mtilde[j]
            for j in keep_fwd
            if sum(i) + sum(j) <= p
        )

    L_rebuilt = {
        n: sum(coeff * Ltilde[i] for i, coeff in terms.items())
        for n, terms in dec.items()
    }

    scale = max(abs(v) for v in L.values())
    for n in full:
        assert abs(L[n] - L_rebuilt[n]) <= 1e-13 * scale * max(1.0, 2 ** sum(n)), (
            p,
            n,
            L[n],
            L_rebuilt[n],
        )
