"""
Stage 2: planar operators against the general (uncompressed, non-planar)
path, evaluated at z=0.

Each planar operator must be numerically equivalent to the general one when
the general one is fed in-plane coordinates (z=0) and, for the multipole
array, values consistent with a source confined to that plane. Unlike the
harmonic-compressed variants (see test_harmonic_operators.py), excluded
multipole entries are provably zero rather than a linear combination of
retained ones, and excluded local entries are simply never read -- so the
translation between compressed and full representations only ever needs two
mechanical primitives (lift-with-zeros, restrict), never the harmonic
relation's Dec-weighted project/expand.

Pure Python; no compiler needed.
"""
import numpy as np
import pytest
import sympy as sp

from fmmgen import generator as g
from fmmgen.utils import Nterms, generate_mappings

symbols = sp.symbols("x y z")
x, y, z = symbols

ORDERS = [1, 2, 3, 4, 5, 6]
TOL = 1e-9

# source_order 0 (charges), 1 (point dipoles, the micromagnetics case) and 2
# (quadrupoles) -- unlike compress, planar imposes no source_order cap, since
# a source's own moment shell (up to n_z = source_order) is exactly what the
# retained set is built to include, not something that itself compresses.
CASES = [(p, s) for s in (0, 1, 2) for p in ORDERS if p > s]


def _setup(p, s=0):
    M_dict, _ = generate_mappings(p, symbols, "grevlex", source_order=s)
    L_dict, _ = generate_mappings(p - s, symbols, "grevlex", source_order=0)
    keep_M, keep_L = g._keep_sets_planar(p, symbols, s, field=True)
    return M_dict, L_dict, keep_M, keep_L


def _inplane_coords(rng):
    """Random x, y with z pinned to 0 -- the defining planar constraint.

    Includes Rinv = 1/R for expressions built from Phi_derivatives (M2L's
    derivs, in particular), which emits the atomic symbol Rinv rather than
    1/sqrt(...); harmless as an unused substitution elsewhere.
    """
    vx, vy = rng.standard_normal(2) * 0.3
    # Keep the point away from the origin: R=0 makes Rinv undefined.
    vx += 0.5 if vx >= 0 else -0.5
    R = (vx**2 + vy**2) ** 0.5
    return {x: float(vx), y: float(vy), z: 0.0, sp.Symbol("Rinv"): 1.0 / R}


def _evaluate(exprs, subs):
    return np.array([float(sp.N(sp.sympify(e).xreplace(subs))) for e in exprs])


def _restrict(vals, full, keep):
    return np.array([vals[full[i]] for i in keep])


def _lift_subs(name, p_or_size, full, keep, comp_vals, is_size=False):
    """Full-indexed substitution putting compressed values in retained slots,
    zero elsewhere -- the planar analogue of harmonic's lift_multipole_subs,
    with the eliminated entries forced to zero rather than left free, since
    that is the whole content of the planar claim being tested."""
    size = p_or_size if is_size else Nterms(p_or_size)
    F = sp.MatrixSymbol(name, size, 1)
    return {
        F[idx]: (float(comp_vals[keep[n]]) if n in keep else 0.0)
        for n, idx in full.items()
    }


def _comp_subs(name, comp_vals):
    C = sp.MatrixSymbol(name, len(comp_vals), 1)
    return {C[i]: float(v) for i, v in enumerate(comp_vals)}


# --------------------------------------------------------------------------

@pytest.mark.parametrize("p,s", CASES)
def test_s2m(p, s):
    """Planar S2M == general S2M restricted to n_z <= s, at z=0."""
    rng = np.random.default_rng(100 + p + 50 * s)
    full, L_dict, keep_M, keep_L = _setup(p, s)
    coords = _inplane_coords(rng)

    nsrc = len([n for n in full if sum(n) == s])
    S = sp.MatrixSymbol("S", nsrc, 1)
    s_subs = {S[i]: float(v) for i, v in enumerate(rng.standard_normal(nsrc))}

    ref = _evaluate(g.generate_S2M_operators(p, symbols, full, source_order=s),
                    {**coords, **s_subs})
    got = _evaluate(
        g.generate_S2M_operators_planar(p, symbols, full, keep_M, source_order=s),
        {**coords, **s_subs},
    )
    assert np.allclose(got, _restrict(ref, full, keep_M), atol=TOL)

    # And the excluded entries of the general operator, evaluated at z=0,
    # really are zero -- not just unused. This is the actual planar claim;
    # the size/index bookkeeping above would pass even if it were false.
    excluded = [ref[full[n]] for n in full if n not in keep_M]
    assert np.allclose(excluded, 0.0, atol=TOL)


@pytest.mark.parametrize("p,s", CASES)
def test_m2m(p, s):
    """Planar M2M agrees after lifting the input and restricting the output."""
    rng = np.random.default_rng(200 + p + 50 * s)
    full, L_dict, keep_M, keep_L = _setup(p, s)
    coords = _inplane_coords(rng)
    Mt = rng.standard_normal(len(keep_M))

    ref = _evaluate(
        g.generate_M_shift_operators(p, symbols, full, source_order=s),
        {**coords, **_lift_subs("M", p, full, keep_M, Mt)},
    )
    got = _evaluate(
        g.generate_M_shift_operators_planar(p, symbols, full, keep_M, source_order=s),
        {**coords, **_comp_subs("M", Mt)},
    )
    assert np.allclose(got, _restrict(ref, full, keep_M), atol=TOL)


@pytest.mark.parametrize("p,s", CASES)
def test_m2l(p, s):
    """Planar M2L (derivs + L) agrees with the general operator at z=0."""
    rng = np.random.default_rng(300 + p + 50 * s)
    full, L_dict, keep_M, keep_L = _setup(p, s)
    coords = _inplane_coords(rng)
    Mt = rng.standard_normal(len(keep_M))

    derivs_full = g.generate_derivs(p, symbols, full, s, harmonic_derivs=False)
    # expansions.L always declares D as Nterms(order)-sized (see expansions.py
    # ~line 169), indexed by M_dict[...] -- a dense but SUBSET range of that
    # declared size whenever s > 0, since M_dict itself is densely reindexed
    # from 0 (generate_mappings). generate_derivs' loop is over the same
    # M_dict, in the same order, so enumerate(derivs_full) lines up with
    # M_dict's values directly.
    D_full = sp.MatrixSymbol("D", Nterms(p), 1)
    Dref = _evaluate(derivs_full, coords)
    d_subs = {D_full[i]: float(v) for i, v in enumerate(Dref)}

    ref = _evaluate(
        g.generate_L_operators(p, symbols, full, L_dict, source_order=s),
        {**d_subs, **_lift_subs("M", p, full, keep_M, Mt)},
    )

    derivs_planar, L_planar = g.generate_M2L_planar(
        p, symbols, full, L_dict, keep_M, keep_L, source_order=s, harmonic_derivs=False,
    )
    Dgot = _evaluate(derivs_planar, coords)
    # Same Nterms(p)-sized D as D_full above: generate_M2L_planar's L_ops
    # comes from the same (unmodified) generate_L_operators/expansions.L, so
    # its embedded D[...] references have the identical declared shape.
    D_planar_sym = sp.MatrixSymbol("D", Nterms(p), 1)
    dp_subs = {D_planar_sym[i]: float(v) for i, v in enumerate(Dgot)}
    got = _evaluate(L_planar, {**dp_subs, **_comp_subs("M", Mt)})

    assert np.allclose(got, _restrict(ref, L_dict, keep_L), atol=TOL)


@pytest.mark.parametrize("p,s", CASES)
def test_l2l(p, s):
    """Planar L2L agrees after lifting the input and restricting the output."""
    rng = np.random.default_rng(400 + p + 50 * s)
    full, L_dict, keep_M, keep_L = _setup(p, s)
    coords = _inplane_coords(rng)
    Lt = rng.standard_normal(len(keep_L))

    ref = _evaluate(
        g.generate_L_shift_operators(p - s, symbols, L_dict, source_order=0),
        {**coords, **_lift_subs("L", p - s, L_dict, keep_L, Lt)},
    )
    got = _evaluate(
        g.generate_L_shift_operators_planar(p - s, symbols, L_dict, keep_L, source_order=0),
        {**coords, **_comp_subs("L", Lt)},
    )
    assert np.allclose(got, _restrict(ref, L_dict, keep_L), atol=TOL)


@pytest.mark.parametrize("p,s", CASES)
def test_l2p(p, s):
    """Potential and field from a planar local expansion match the general one
    evaluated at a target also in the z=0 plane."""
    rng = np.random.default_rng(500 + p + 50 * s)
    full, L_dict, keep_M, keep_L = _setup(p, s)
    coords = _inplane_coords(rng)
    Lt = rng.standard_normal(len(keep_L))
    q = p - s

    ref = _evaluate(
        g.generate_L2P_operators(q, symbols, L_dict),
        {**coords, **_lift_subs("L", q, L_dict, keep_L, Lt)},
    )
    got = _evaluate(
        g.generate_L2P_operators_planar(q, symbols, L_dict, keep_L),
        {**coords, **_comp_subs("L", Lt)},
    )
    assert np.allclose(got, ref, atol=TOL)
    # Fz (component index 3, when field is on) is exactly zero for a target
    # also confined to z=0: the potential's only z-content is odd powers of
    # the L array's own n_z index, and evaluating a n_z<=1 local expansion's
    # z-derivative at dz=0 kills the n_z=0 shell while the n_z=1 shell's own
    # coefficient contributes a genuinely nonzero constant term -- so this is
    # NOT asserted to be zero here (Lt[keep_L for n_z=1] is random and
    # nonzero); the physical zero-Fz case (monopole source, both source and
    # target in-plane) is instead demonstrated end-to-end in test_planar.py.


@pytest.mark.parametrize("p,s", CASES)
def test_m2p(p, s):
    """Potential and field from a planar multipole match the general one."""
    rng = np.random.default_rng(600 + p + 50 * s)
    full, L_dict, keep_M, keep_L = _setup(p, s)
    # M2P differentiates 1/R inline, so keep the target well separated from
    # the origin, and supply Rinv consistently -- Phi_derivatives emits the
    # atomic symbol Rinv rather than 1/sqrt(...).
    coords = {x: 1.7, y: -2.3, z: 0.0}
    coords[sp.Symbol("Rinv")] = 1.0 / np.sqrt(1.7**2 + 2.3**2)
    Mt = rng.standard_normal(len(keep_M))

    ref = _evaluate(
        g.generate_M2P_operators(p, symbols, full, source_order=s),
        {**coords, **_lift_subs("M", p, full, keep_M, Mt)},
    )
    got = _evaluate(
        g.generate_M2P_operators_planar(p, symbols, full, keep_M, source_order=s),
        {**coords, **_comp_subs("M", Mt)},
    )
    assert np.allclose(got, ref, atol=TOL)


@pytest.mark.parametrize("s", (0, 1, 2))
def test_p2p(s):
    """
    Planar P2P matches the general one at z=0, for both particles of the
    pair. Unlike every other planar operator there is no array to
    lift/restrict here -- P2P reads the caller's raw source array directly,
    sized by source_order alone -- so this exercises generate_P2P_operators_planar
    on its own rather than through _setup/_keep_sets_planar.
    """
    rng = np.random.default_rng(700 + 50 * s)
    M_dict, _ = generate_mappings(s, symbols, source_order=s)
    ops = g.generate_P2P_operators(symbols, M_dict, potential=True, field=True, source_order=s)
    ops_xy = g.generate_P2P_operators_planar(symbols, M_dict, potential=True, field=True,
                                             source_order=s, ops=ops)
    assert len(ops) == len(ops_xy)  # output size is NOT reduced -- see the
    # module docstring on generate_P2P_operators_planar for why Fz is kept.

    nsrc = Nterms(s)  # generate_P2P_operators' own internal S array size
    S = sp.MatrixSymbol("S", nsrc, 1)
    s_subs = {S[i]: float(v) for i, v in enumerate(rng.standard_normal(nsrc))}
    dx, dy = float(rng.uniform(0.5, 2)), float(rng.uniform(0.5, 2))
    coords = {x: dx, y: dy, z: 0.0, sp.Symbol("Rinv"): 1.0 / np.sqrt(dx**2 + dy**2)}

    ref = _evaluate(ops, {**coords, **s_subs})
    got = _evaluate(ops_xy, {**coords, **s_subs})
    assert np.allclose(got, ref, atol=TOL)

    if s == 0:
        # The one case where Fz (last entry, potential+field both on) really
        # is exactly zero: a monopole's potential is even in z, so its
        # z-derivative vanishes at z=0. Pinned symbolically, not just at one
        # random point -- see test_planar.py for the s=1 contrast (genuinely
        # nonzero, via the dipole's own muz).
        assert sp.simplify(sp.sympify(ops_xy[-1])) == 0
