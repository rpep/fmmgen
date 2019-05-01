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
import functools


q, x, y, z = sp.symbols('q x y z')

def TriangleNumbers(n):
    """
    Returns the nth triangle number
    """
    return int(n * (n + 1) / 2)

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
    monom_key = monomial_key('grevlex', symbols)
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




"""Tools and arithmetics for monomials of distributed polynomials. """

from itertools import combinations_with_replacement, product
from sympy.core import Mul, S

#
#
# def itermonomials(variables, degree):
#     """
#     Generate a set of monomials of the given total degree or less.
#
#     Given a set of variables `V` and a total degree `N` generate
#     a set of monomials of degree at most `N`. The total number of
#     monomials in commutative variables is huge and is given by the
#     following formula:
#
#     .. math::
#
#         \frac{(\#V + N)!}{\#V! N!}
#
#     For example if we would like to generate a dense polynomial of
#     a total degree `N = 50` in 5 variables, assuming that exponents
#     and all of coefficients are 32-bit long and stored in an array we
#     would need almost 80 GiB of memory! Fortunately most polynomials,
#     that we will encounter, are sparse.
#
#     Examples
#     ========
#
#     Consider monomials in commutative variables `x` and `y`
#     and non-commutative variables `a` and `b`::
#
#         >>> from sympy import symbols
#         >>> from sympy.polys.monomials import itermonomials
#         >>> from sympy.polys.orderings import monomial_key
#         >>> from sympy.abc import x, y
#
#         >>> sorted(itermonomials([x, y], 2), key=monomial_key('grlex', [y, x]))
#         [1, x, y, x**2, x*y, y**2]
#
#         >>> sorted(itermonomials([x, y], 3), key=monomial_key('grlex', [y, x]))
#         [1, x, y, x**2, x*y, y**2, x**3, x**2*y, x*y**2, y**3]
#
#         >>> a, b = symbols('a, b', commutative=False)
#         >>> itermonomials([a, b, x], 2)
#         {1, a, a**2, b, b**2, x, x**2, a*b, b*a, x*a, x*b}
#
#
#     """
#     if degree < 0:
#         return set()
#     if not variables or degree == 0:
#         return [S(1)]
#     # Force to list in case of passed tuple or other incompatible collection
#     variables = list(variables) + [S(1)]
#     if all(variable.is_commutative for variable in variables):
#         return [Mul(*item) for item in combinations_with_replacement(variables, degree)]
#     else:
#         return [Mul(*item) for item in product(variables, repeat=degree)]
