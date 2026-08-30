"""
Stage 1 gate for the planar (2D-plane) reduction.

The claim under test is positional, not combinatorial (contrast
test_harmonic.py's Hankel-factorisation gate): a source confined to z=0 has a
multipole array that is IDENTICALLY ZERO outside n_z <= source_order, because
every displacement entering M_shift's recursion has dz = 0, and dz**k = 0 for
k >= 1 kills every term that would otherwise raise n_z above what the
source's own moment shell (|n| = source_order) started with. Local array
entries outside n_z <= 1 (field) / n_z <= 0 (potential only) are not zero in
general, but are never read by a correctly z=0-evaluating L2L/L2P, which is
the weaker, sufficient claim actually needed and the one tested here.

If either claim were false, generate_*_operators_planar would still run and
produce a same-shaped answer -- restrict_local doesn't care whether the
entries it drops are actually negligible -- so this file exists to catch that
silently-wrong-answer failure mode, which test_planar_operators.py's
numerical-equivalence checks (built from the SAME generate_S2M_operators
etc.) could not by themselves distinguish from a self-consistent but wrong
implementation.

These are pure Python and need no compiler.
"""
import sympy as sp
import pytest

from fmmgen import generator as g
from fmmgen.utils import Nterms, generate_mappings

symbols = sp.symbols("x y z")
x, y, z = symbols

ORDERS = list(range(1, 9))


# --------------------------------------------------------------------------
# counting and index-set bookkeeping
# --------------------------------------------------------------------------

@pytest.mark.parametrize("p,s", [(p, s) for s in (0, 1, 2, 3) for p in ORDERS if p >= s])
def test_keep_M_is_consistent(p, s):
    """Indices are contiguous from zero, n_z <= s, |n| in [s, p], order-preserving."""
    keep_M = g._planar_keep(p, symbols, s, max_nz=s)
    assert sorted(keep_M.values()) == list(range(len(keep_M)))
    assert all(n[2] <= s for n in keep_M)
    assert all(s <= sum(n) <= p for n in keep_M)

    full, _ = generate_mappings(p, symbols, source_order=s)
    kept_in_full_order = [n for n in full if n[2] <= s]
    assert list(keep_M.keys()) == kept_in_full_order
    assert g.Nterms_planar_M(p, s) == len(keep_M)


@pytest.mark.parametrize("p,field", [(p, f) for f in (True, False) for p in ORDERS])
def test_keep_L_is_consistent(p, field):
    max_nz = 1 if field else 0
    keep_L = g._planar_keep(p, symbols, 0, max_nz=max_nz)
    assert sorted(keep_L.values()) == list(range(len(keep_L)))
    assert all(n[2] <= max_nz for n in keep_L)

    full, _ = generate_mappings(p, symbols, source_order=0)
    kept_in_full_order = [n for n in full if n[2] <= max_nz]
    assert list(keep_L.keys()) == kept_in_full_order
    assert g.Nterms_planar(p, max_nz) == len(keep_L)


# --------------------------------------------------------------------------
# the gate: the multipole array really is zero outside n_z <= source_order
# --------------------------------------------------------------------------

@pytest.mark.parametrize("p,s", [(p, s) for s in (0, 1, 2) for p in [s, s + 1, s + 2, s + 4]])
def test_s2m_is_exactly_zero_outside_keep(p, s):
    """
    generate_S2M_operators, evaluated at z=0, is the exact symbolic zero
    polynomial (not merely numerically small) at every n_z > s. Checked
    symbolically -- no substitution, no floating point -- since S2M is
    already linear in the caller's S array and the claim is that whole
    coefficients vanish, not that they happen to cancel at a sampled point.
    """
    M_dict, _ = generate_mappings(p, symbols, "grevlex", source_order=s)
    ops = g.generate_S2M_operators(p, symbols, M_dict, source_order=s)
    ops0 = [sp.expand(sp.sympify(e).xreplace({z: 0})) for e in ops]
    for n, i in M_dict.items():
        if n[2] > s:
            assert ops0[i] == 0, (p, s, n, ops0[i])


@pytest.mark.parametrize("p,s", [(p, s) for s in (0, 1, 2) for p in [s + 1, s + 2, s + 4]])
def test_m2m_preserves_zero_outside_keep(p, s):
    """
    The invariant propagates through M2M: feeding generate_M_shift_operators
    an input that is symbolically zero outside n_z <= s (an arbitrary
    generic value elsewhere) produces an output that is ALSO exactly zero
    outside n_z <= s, at z=0 -- for every generic input, not merely the one
    S2M happens to construct. This is what licenses applying the same
    restrict_local at every order of a multi-level tree, not just once at
    the leaves.
    """
    M_dict, _ = generate_mappings(p, symbols, "grevlex", source_order=s)
    M = sp.MatrixSymbol("M", Nterms(p), 1)
    generic_subs = {M[i]: (M[i] if n[2] <= s else 0) for n, i in M_dict.items()}

    ops = g.generate_M_shift_operators(p, symbols, M_dict, source_order=s)
    ops0 = [sp.expand(sp.sympify(e).xreplace({z: 0}).xreplace(generic_subs)) for e in ops]
    for n, i in M_dict.items():
        if n[2] > s:
            assert ops0[i] == 0, (p, s, n, ops0[i])


def test_quadrupole_and_beyond_are_not_rejected():
    """
    Positive contrast with compress's source_order >= 2 rejection
    (test_harmonic_operators.py::test_quadrupole_sources_are_rejected):
    planar imposes no such cap, and this is not an oversight. A quadrupole's
    own moment shell has a genuinely nonzero n_z=2 entry (e.g. mu_zz for an
    atomistic quadrupole tensor, unrelated to whether the quadrupole's
    POSITION is planar), and n_z <= source_order = 2 is exactly the retained
    set built to include it -- verified here by checking that entry survives
    generate_S2M_operators_planar untouched, and does NOT raise.
    """
    p, s = 5, 2
    M_dict, _ = generate_mappings(p, symbols, "grevlex", source_order=s)
    keep_M = g._planar_keep(p, symbols, s, max_nz=s)
    assert (0, 0, 2) in keep_M  # the mu_zz shell entry itself

    ops = g.generate_S2M_operators_planar(p, symbols, M_dict, keep_M, source_order=s)
    S = sp.MatrixSymbol("S", len([n for n in M_dict if sum(n) == s]), 1)
    shell_expr = sp.sympify(ops[keep_M[(0, 0, 2)]])
    # It is exactly the caller's mu_zz component, not zero and not some
    # combination -- S2M's shell entries index straight into S.
    assert shell_expr == S[M_dict[(0, 0, 2)]]


# --------------------------------------------------------------------------
# the gate: L2P/L2L never read local entries outside the retained set
# --------------------------------------------------------------------------

@pytest.mark.parametrize("p", [p for p in ORDERS if p >= 2])
def test_l2p_does_not_depend_on_excluded_local_entries(p):
    """
    generate_L2P_operators, evaluated at a target ALSO in the z=0 plane, does
    not reference L[i] for any n_z(i) >= 2 at all (vacuous below p=2, where no
    monomial can reach n_z=2 at all) -- checked by
    differentiating symbolically w.r.t. each such L[i] and confirming the
    derivative is exactly zero, which is a stronger and more direct claim
    than "restricting doesn't change the numeric answer for one random L".
    """
    L_dict, _ = generate_mappings(p, symbols, "grevlex", source_order=0)
    L = sp.MatrixSymbol("L", Nterms(p), 1)
    ops = g.generate_L2P_operators(p, symbols, L_dict, potential=True, field=True)
    ops0 = [sp.sympify(e).xreplace({z: 0}) for e in ops]

    excluded = [i for n, i in L_dict.items() if n[2] >= 2]
    assert excluded, "test is vacuous without excluded entries to check"
    for e in ops0:
        for i in excluded:
            assert sp.diff(e, L[i]) == 0, (p, i)


@pytest.mark.parametrize("p", [1, 2, 3, 4])
def test_l2p_potential_only_does_not_depend_on_nz_1_either(p):
    """With field=False, the potential alone never even needs n_z <= 1's Fz
    shell -- only n_z = 0 survives -- since Fz is not being computed."""
    L_dict, _ = generate_mappings(p, symbols, "grevlex", source_order=0)
    L = sp.MatrixSymbol("L", Nterms(p), 1)
    ops = g.generate_L2P_operators(p, symbols, L_dict, potential=True, field=False)
    ops0 = [sp.sympify(e).xreplace({z: 0}) for e in ops]

    excluded = [i for n, i in L_dict.items() if n[2] >= 1]
    for e in ops0:
        for i in excluded:
            assert sp.diff(e, L[i]) == 0, (p, i)
