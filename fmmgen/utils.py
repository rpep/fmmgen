#########################################
#
# Utility functions for the FMM
#
# (C) Ryan Pepper, 2018
# University of Southampton, UK
#
#
#########################################
import sympy as sp
from sympy.polys.orderings import monomial_key
from sympy.polys.monomials import itermonomials as sp_itermonomials
import functools

q, x, y, z = sp.symbols("q x y z")


def TriangleNumbers(n):
    """
    Returns the nth triangle number
    """
    return int((n * (n + 1)) / 2)


@functools.lru_cache(maxsize=None)
def Nterms(p):
    """
    Determines the number of terms in a multipole expansion of order
    n
    """
    if p < 0:
        return 0
    else:
        return int(sum([TriangleNumbers(i) for i in range(p + 2)]))


def new_itermonomials(symbols, lower_order, upper_order):
    monom_key = monomial_key("grevlex", symbols)
    monoms = itermonomials(symbols, upper_order)
    monoms = sorted(monoms, key=monom_key)
    new_monoms = []
    for monom in monoms:
        monom_dict = monom.as_powers_dict()
        order = 0
        for symbol in symbols:
            order += monom_dict[symbol]
        if order >= lower_order:
            new_monoms.append(monom)
    return set(new_monoms)


def itermonomials(symbols, max_degree, min_degree=0):
    monoms = list(sp_itermonomials(symbols, max_degree))
    new_monoms = []
    for monom in monoms:
        monom_dict = monom.as_powers_dict()
        order = 0
        for symbol in symbols:
            order += monom_dict[symbol]
            if order >= min_degree:
                new_monoms.append(monom)
    return set(new_monoms)


def generate_mappings(order, symbols, key="grevlex", source_order=0):
    """
    generate_mappings(order, symbols, key='grevlex'):

    Generates a set of mappings between three-tuples of
    indices denoting monomials and and array indices.
    Returns both the forward and backward maps.

    Inputs:
    order: int
        Maximum monomial order

    symbols: list
        List of sympy.Symbol type objects

    source_order: int
        Integer describing order of o

    Returns:
    dict:
        Forward mapping from n-tuple to array index.

    dict:
        Reversed version; mapping from array index to
        tuple mapping.

    Example:
    >>> x, y, z = sp.symbols('x y z')
    >>> map, rmap = generate_mappings(1, [x, y, z])
    >>> print(map):
    {(0, 0, 0): 0, (1, 0, 0): 1, (0, 1, 0): 2, (0, 0, 1): 3}
    >>> print(rmap):
    {0: (0, 0, 0), 1: (1, 0, 0), 2: (0, 1, 0), 3: (0, 0, 1)}
    """
    # Cached on (order, symbols, key, source_order): every generator.py
    # operator (M_shift, L, L_shift, ...) calls this once per multi-index n
    # it emits, with the SAME order/symbols/source_order each time, so
    # without caching the itermonomials + grevlex sort was being rebuilt from
    # scratch O(Nterms(order)) times per operator. Safe because no caller
    # mutates the returned dicts (verified: every dict subscript-assignment
    # elsewhere in the package targets a locally-built dict, never the
    # index_dict/M_dict/L_dict returned from here).
    return _generate_mappings_cached(order, tuple(symbols), key, source_order)


@functools.lru_cache(maxsize=None)
def _generate_mappings_cached(order, symbols, key, source_order):
    if order < source_order:
        raise ValueError("source_order must be <= order for meaningful calculations to occur")

    x, y, z = symbols
    rsymbols = [z, y, x]

    monoms = itermonomials(symbols, order, source_order)
    if key:
        monom_key = monomial_key(key, rsymbols)
        monoms = sorted(monoms, key=monom_key)

    index_dict = {}
    rindex_dict = {}
    for i, monom in enumerate(monoms):
        d = monom.as_powers_dict()
        n = d[x], d[y], d[z]
        index_dict[n] = i
        rindex_dict[i] = n
    return index_dict, rindex_dict
