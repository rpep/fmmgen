import sympy as sp
from .expansions import M, M_shift, L, L_shift, phi_deriv, Phi_derivatives
from .utils import generate_mappings, Nterms
import logging

logger = logging.getLogger(name="fmmgen")


def generate_M_operators(order, symbols, M_dict):
    """
    generate_M_operators(order, symbols, index_dict):

    Generates multipole operators up to order.

    Input:
    order, int:
        Maximum order of multipole expansion

    symbols, list:
        List of sympy symbol type objects which
        define coordinate labels.

    index_dict:
        Forward mapping dictionary between
        monomials of symbols and array indices,
        generated by generate_mappings or otherwise.

    Output:
    list:
        List of symbolic multipole moments up to order.

    Example:
    >>> order = 2
    >>> x, y, z = sp.symbols('x y z')
    >>> map, _ = generate_mappings(order, [x, y, z])
    >>> generate_M_operators(order, (x, y, z), map)
    [q, -q*x, -q*y, -q*z, q*x**2/2, q*x*y, q*x*z, q*y**2/2, q*y*z, q*z**2/2]
    """
    x, y, z = symbols
    M_operators = []
    for n in M_dict.keys():
        M_operators.append(M(n, symbols))
    return M_operators


def generate_M_shift_operators(order, symbols, M_dict, source_order=0):
    """
    generate_M_shift_operators(order, symbols, index_dict):

    Generates multipole shifting operators up to order.

    Input:
    order, int:
        Maximum order of multipole expansion

    symbols, list:
        List of sympy symbol type objects which define coordinate labels.

    index_dict:
        Forward mapping dictionary between
        monomials of symbols and array indices, generated by generate_mappings
         or otherwise.

    Output:
    list:
        List of symbolic multipole shifting operators up to order.

    Example:
    >>> order = 1
    >>> x, y, z = sp.symbols('x y z')
    >>> map, _ = generate_mappings(order, [x, y, z])
    >>> generate_M_shift_operators(order, (x, y, z), map)
    [M[0, 0], x*M[0, 0] + M[1, 0], y*M[0, 0] + M[2, 0], z*M[0, 0] + M[3, 0]]
    """
    x, y, z = symbols
    M_operators = []
    for n in M_dict.keys():
        M_operators.append(M_shift(n, order, symbols, M_dict, source_order=source_order))
    return M_operators


def generate_derivs(order, symbols, M_dict, source_order=0, harmonic_derivs=False):
    D = sp.MatrixSymbol("D", Nterms(order), 1)
    derivs = []
    for n in M_dict.keys():
        if n[2] > 1 and harmonic_derivs:
            k = (n[0], n[1], n[2] - 2)
            k1 = (k[0] + 2, k[1], k[2])
            k2 = (k[0], k[1] + 2, k[2])
            derivs.append(-D[M_dict[k1]] - D[M_dict[k2]])
        else:
            derivs.append(Phi_derivatives(n, symbols))
    return derivs


def generate_L_operators(order, symbols, M_dict, L_dict, source_order=0):
    """
    generate_L_operators(order, symbols, index_dict):

    Generates local expansion operators up to given order.

    Input:
    order, int:
        Maximum order of multipole expansion

    symbols, list:
        List of sympy symbol type objects which define coordinate labels.

    index_dict:
        Forward mapping dictionary between
        monomials of symbols and array indices, generated by generate_mappings
        or otherwise.

    inline_derivs: bool
        Calculate derivatives inline rather than precalculating these.

    Output:
    list:
        List of symbolic local expansion operators up to order.

    Example:
    >>> order = 1
    >>> x, y, z = sp.symbols('x y z')
    >>> map, _ = generate_mappings(order, (x, y, z))
    >>> generate_L_operators(order, (x, y, z), map)
    [M[0, 0]/R - 1.0*x*M[1, 0]/R**3 - 1.0*y*M[2, 0]/R**3 - 1.0*z*M[3, 0]/R**3,
    -1.0*x*M[0, 0]/R**3, -1.0*y*M[0, 0]/R**3, -1.0*z*M[0, 0]/R**3]
    """
    x, y, z = symbols
    L_operators = []
    for n in L_dict.keys():
        L_operators.append(L(n, order, symbols, M_dict, source_order=source_order, eval_derivs=False))

    return L_operators


def generate_L_shift_operators(order, symbols, L_dict, source_order=0):
    """
    generate_L_shift_operators(order, symbols, index_dict):

    Generates multiple operators up to order.

    Input:
    order, int:
        Maximum order of multipole expansion

    symbols, list:
        List of sympy symbol type objects which
        define coordinate labels.

    index_dict:
        Forward mapping dictionary between
        monomials of symbols and array indices,
        generated by generate_mappings or otherwise.

    Output:
    list:
        List of symbolic local expansion shifting operators up to order.

    Example:
    >>> order = 1
    >>> x, y, z = sp.symbols('x y z')
    >>> map, _ = generate_mappings(order, (x, y, z))
    >>> generate_L_shift_operators(order, (x, y, z), map)
    [x*L[1, 0] + y*L[2, 0] + z*L[3, 0] + L[0, 0], L[1, 0], L[2, 0], L[3, 0]]
    """
    x, y, z = symbols
    L_shift_operators = []
    for n in L_dict.keys():
        L_shift_operators.append(L_shift(n, order, symbols, L_dict, source_order=source_order))
    return L_shift_operators


def generate_M2P_operators(
    order,
    symbols,
    M_dict,
    potential=True,
    field=True,
    source_order=0,
    harmonic_derivs=False,
):
    """
    generate_M2L_operators(order, symbols, index_dict)

    Generates potential and field calculation operators for the
    Barnes-Hut method up to order.
    """
    x, y, z = symbols
    R = (x**2 + y**2 + z**2) ** 0.5

    terms = []

    V = L((0, 0, 0), order, symbols, M_dict, source_order=source_order, eval_derivs=True).subs("R", R)
    if potential:
        terms.append(V.subs(1 / R, "Rinv"))

    if field:
        Fx = -sp.diff(V, x).subs(1 / R, "Rinv")
        Fy = -sp.diff(V, y).subs(1 / R, "Rinv")
        Fz = -sp.diff(V, z).subs(1 / R, "Rinv")
        terms.append(Fx)
        terms.append(Fy)
        terms.append(Fz)

    return terms


def generate_L2P_operators(order, symbols, L_dict, potential=True, field=True):
    """
    generate_L2P_operators(order, symbols, index_dict):

    Generates potential and field calculation operators for the Fast
    Multipole Method up to order.

    Input:
    order, int:
        Maximum order of multipole expansion

    symbols, list:
        List of sympy symbol type objects which define coordinate labels.

    index_dict:
        Forward mapping dictionary between monomials of symbols and array
        indices, generated by generate_mappings or otherwise.

    Output:
    list:
        List of symbolic field calculation operators from local expansions up to
        order.

    Example:
    >>> order = 1
    >>> x, y, z = sp.symbols('x y z')
    >>> map, _ = generate_mappings(order, (x, y, z))
    >>> generate_L2P_operators(order, (x, y, z), map)
    [x*L[1, 0] + y*L[2, 0] + z*L[3, 0] + L[0, 0], -L[1, 0], -L[2, 0], -L[3, 0]]
    """
    x, y, z = symbols

    terms = []

    if potential:
        V = phi_deriv(order, symbols, L_dict, deriv=(0, 0, 0))
        terms.append(V)

    if field:
        Fx = -phi_deriv(order, symbols, L_dict, deriv=(1, 0, 0))
        Fy = -phi_deriv(order, symbols, L_dict, deriv=(0, 1, 0))
        Fz = -phi_deriv(order, symbols, L_dict, deriv=(0, 0, 1))
        terms.append(Fx)
        terms.append(Fy)
        terms.append(Fz)
    return terms


def generate_P2P_operators(symbols, M_dict, potential=True, field=True, source_order=0):
    order = source_order
    M_dict, _ = generate_mappings(source_order, symbols, "grevlex", source_order=source_order)
    x, y, z = sp.symbols("x y z")
    R = (x**2 + y**2 + z**2) ** 0.5

    S_map, _ = generate_mappings(source_order, [x, y, z], key="grevlex", source_order=source_order)
    # print('S_map = {}'.format(S_map))

    M = sp.MatrixSymbol("M", Nterms(order), 1)
    S = sp.MatrixSymbol("S", Nterms(order), 1)

    subsdict = {M[i]: 0 for i in range(Nterms(order))}

    for key in S_map.keys():
        subsdict[M[M_dict[key]]] = S[M_dict[key]]

    V = L((0, 0, 0), order, symbols, M_dict, source_order=source_order).subs("R", R).subs(subsdict)

    terms = []
    # Note: R must be substituted late for correct derivatives!
    if potential:
        terms.append(V.subs(1 / R, "Rinv"))

    if field:
        Fx = -sp.diff(V, x).subs(1 / R, "Rinv")
        Fy = -sp.diff(V, y).subs(1 / R, "Rinv")
        Fz = -sp.diff(V, z).subs(1 / R, "Rinv")
        terms.append(Fx)
        terms.append(Fy)
        terms.append(Fz)
    return terms
