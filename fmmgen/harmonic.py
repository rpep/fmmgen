#########################################
#
# Harmonic (traceless) compression of Cartesian Taylor expansions
#
#########################################
r"""
Compression of Cartesian Taylor expansions using harmonicity of 1/R.

The Laplace Green's function is harmonic away from the origin, so its
derivatives are trace-free:

    \nabla^2 D_k = D_{k+(2,0,0)} + D_{k+(0,2,0)} + D_{k+(0,0,2)} = 0.

Under fmmgen's grevlex ordering on [z, y, x] the natural pivot is the z index,
so every index with n_z >= 2 can be eliminated by recursive application of

    D_{k+(0,0,2)} = -D_{k+(2,0,0)} - D_{k+(0,2,0)}                       (*)

leaving the retained set

    keep(p) = {n : |n| <= p, n_z <= 1}

of size C(p+2,2) + C(p+1,2) = (p+1)^2, against Nterms(p) = C(p+3,3) for the
uncompressed array -- a factor p/3 fewer coefficients, 2.53x at p = 11.

Two properties of (*) matter and are exploited elsewhere:

1. It preserves total degree. k+(2,0,0), k+(0,2,0) and k+(0,0,2) all have
   degree |k|+2, so an elimination never references an index of higher degree
   than the one it replaces. This is what makes the compression compatible
   with the |n| + |m| <= p truncation used throughout fmmgen.

2. Repeated application collects rather than branching: eliminating n_z = 2m
   yields m+1 distinct terms, not 2^m. For example
   D_{(0,0,4)} = D_{(4,0,0)} + 2 D_{(2,2,0)} + D_{(0,4,0)}, which is just
   \partial_z^4 = (\partial_x^2 + \partial_y^2)^2.

This is the Laplace instance of the PDE-driven compression of Fernando and
Klockner (arXiv:2305.17867, 2023). The identity is much older: Applequist
(1989), Hinsen and Felderhof (1992); Dehnen (2014, Appendix A.2) names it in
an FMM context. fmmgen has applied (*) to the derivative array since the
harmonic_derivs option was added -- see generate_derivs in generator.py -- but
not to the stored multipole and local coefficients, which is where the
asymptotic saving lives.

Only the Laplace kernel is handled, because expansions.py hardcodes phi = 1/R.
If that is ever generalised, keep_mappings and decompress are the two functions
that would need to take the PDE as an argument.
"""
import functools

import sympy as sp

from .utils import Nterms, generate_mappings


def Nkeep(p):
    """
    Number of retained coefficients in a harmonic expansion of order p.

    >>> [Nkeep(p) for p in range(5)]
    [1, 4, 9, 16, 25]
    """
    if p < 0:
        return 0
    return (p + 1) ** 2


def keep_mappings(order, symbols, key="grevlex", source_order=0):
    """
    keep_mappings(order, symbols, key='grevlex', source_order=0)

    As utils.generate_mappings, but restricted to the retained multi-indices
    {n : |n| <= order, n_z <= 1}. Ordering is inherited from
    generate_mappings, so retained indices appear in the same relative order
    as in the uncompressed array and degree shells remain contiguous.

    Returns the same (forward, reverse) pair of dicts, re-indexed to be
    contiguous from zero.

    Example:
    >>> import sympy as sp
    >>> x, y, z = sp.symbols('x y z')
    >>> fwd, rev = keep_mappings(2, [x, y, z])
    >>> fwd
    {(0, 0, 0): 0, (1, 0, 0): 1, (0, 1, 0): 2, (0, 0, 1): 3, (2, 0, 0): 4, \
(1, 1, 0): 5, (1, 0, 1): 6, (0, 2, 0): 7, (0, 1, 1): 8}
    """
    full, _ = generate_mappings(order, symbols, key=key, source_order=source_order)
    # dicts preserve insertion order, and generate_mappings inserts in sorted
    # monomial order, so filtering keeps the relative ordering intact.
    kept = [n for n in full if n[2] <= 1]
    index_dict = {n: i for i, n in enumerate(kept)}
    rindex_dict = {i: n for i, n in enumerate(kept)}
    return index_dict, rindex_dict


def pivot_step(n):
    """
    pivot_step(n)

    One application of the trace-free relation: the pair (k1, k2) with

        D[n] = -D[k1] - D[k2]

    or None if n is already retained (n_z <= 1).

    k1 and k2 have the same total degree as n but z index reduced by two, so
    under the grevlex ordering on [z, y, x] they always appear *earlier* in
    the array than n. That is what lets generate_derivs emit a single
    subtraction per eliminated entry, referring back to entries already
    computed, rather than expanding the full recursion -- contrast
    expand_index, which is needed when the target basis rather than the array
    order is what matters.

    >>> pivot_step((1, 0, 1)) is None
    True
    >>> pivot_step((0, 0, 2))
    ((2, 0, 0), (0, 2, 0))
    """
    if n[2] <= 1:
        return None
    k = (n[0], n[1], n[2] - 2)
    return (k[0] + 2, k[1], k[2]), (k[0], k[1] + 2, k[2])


@functools.lru_cache(maxsize=None)
def expand_index(n):
    """
    expand_index(n)

    Express the multi-index n in terms of retained indices, by recursive
    application of the trace-free relation (*).

    Returns a dict {retained multi-index: integer coefficient}. Indices which
    are already retained (n_z <= 1) map to themselves with coefficient 1.

    Example:
    >>> expand_index((0, 0, 2)) == {(2, 0, 0): -1, (0, 2, 0): -1}
    True
    >>> expand_index((0, 0, 4)) == {(4, 0, 0): 1, (2, 2, 0): 2, (0, 4, 0): 1}
    True
    """
    step = pivot_step(n)
    if step is None:
        return {n: 1}

    out = {}
    for child in step:
        for idx, coeff in expand_index(child).items():
            total = out.get(idx, 0) - coeff
            if total:
                out[idx] = total
            else:
                out.pop(idx, None)
    return out


def decompress(order, symbols, key="grevlex", source_order=0):
    """
    decompress(order, symbols, key='grevlex', source_order=0)

    The decompression map Dec, as a dict

        {full multi-index: {retained multi-index: coefficient}}

    covering every multi-index with |n| <= order. Retained indices map to
    themselves. Because (*) preserves total degree, every key in the inner
    dict has the same total degree as the outer key.
    """
    full, _ = generate_mappings(order, symbols, key=key, source_order=source_order)
    return {n: expand_index(n) for n in full}


# --------------------------------------------------------------------------
# Turning uncompressed operators into compressed ones
#
# Every operator is linear in its expansion array, so compression is four
# mechanical transformations rather than six rewritten derivations. Which one
# applies depends on whether the array is a multipole or a local, and whether
# it is an input or an output.
#
#   multipole in   lift_multipole_subs   retained slot, or zero
#   multipole out  project_multipole     Dec^T
#   local in       expand_local_subs     Dec
#   local out      restrict_local        take the retained subset
#
# The multipole cases are asymmetric with the local ones, and the reason is
# worth stating. A multipole array is *not* trace-free -- M_n = sum q x^n / n!
# are moments of a source distribution, unconstrained. What is true is that
# the field it represents,
#
#     phi(y) = sum_n M_n D_n(y - z) = sum_{i in keep} (Dec^T M)_i D_i(y - z),
#
# depends on M only through Dec^T M, because D is trace-free and an inner
# product with a trace-free tensor sees only the trace-free part of its
# co-operand. So any representative with the correct projection may be
# substituted, and zeros on the eliminated entries are the cheapest one --
# note Dec^T applied to that lift returns the compressed array unchanged,
# since Dec[n, i] = delta_{ni} for retained n.
#
# A local array *is* constrained: L_n = sum_m M_m D_{n+m} inherits the
# relation from D, so the eliminated entries are determined and restriction
# loses nothing.
# --------------------------------------------------------------------------

def lift_multipole_subs(name, order, full_dict, keep_dict):
    """
    Substitution sending a full-indexed multipole array onto the compressed
    one: retained entries to their compressed slot, eliminated entries to zero.
    """
    full = sp.MatrixSymbol(name, Nterms(order), 1)
    comp = sp.MatrixSymbol(name, Nkeep(order), 1)
    return {
        full[i]: (comp[keep_dict[n]] if n in keep_dict else sp.Integer(0))
        for n, i in full_dict.items()
    }


def expand_local_subs(name, order, symbols, full_dict, keep_dict):
    """
    Substitution sending a full-indexed local array onto the compressed one,
    L_n = sum_i Dec[n, i] Lt_i.
    """
    full = sp.MatrixSymbol(name, Nterms(order), 1)
    comp = sp.MatrixSymbol(name, Nkeep(order), 1)
    dec = decompress(order, symbols)
    return {
        full[i]: sum(c * comp[keep_dict[k]] for k, c in dec[n].items())
        for n, i in full_dict.items()
    }


def project_multipole(exprs, order, symbols, full_dict, keep_dict):
    """Dec^T applied to a full-indexed expression list, giving a keep-indexed one."""
    dec = decompress(order, symbols)
    acc = {i: sp.Integer(0) for i in keep_dict}
    for n, idx in full_dict.items():
        for k, c in dec[n].items():
            acc[k] += c * exprs[idx]
    return [acc[i] for i in keep_dict]


def restrict_local(exprs, full_dict, keep_dict):
    """The retained subset of a full-indexed expression list."""
    return [exprs[full_dict[i]] for i in keep_dict]
