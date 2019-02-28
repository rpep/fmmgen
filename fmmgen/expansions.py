import sympy as sp
from sympy.polys.orderings import monomial_key
from sympy import factorial
from sympy import itermonomials
from .utils import q, Nterms
import functools



def M(n, symbols):
    """
    M(n, symbols)

    Returns symbolic expression for the n = (nx, ny, nz) multipole expansion of
    a charge.

    i.e if we want the x-component of the dipole moment:

    >>> x, y, z = sp.symbols('x y z')
    >>> M((1, 0, 0), (x, y, z))
    -q*x
    """
    dx, dy, dz = symbols
    return q * (-1)**(n[0] + n[1] + n[2]) * dx**n[0] * dy**n[1] * dz**n[2] / \
        (factorial(n[0]) * factorial(n[1]) * factorial(n[2]))


def M_shift(n, order, symbols, index_dict):
    """
    M_shift(n, order, symbols, index_dict)

    Returns symbolic expression for the n = (nx, ny, nz) multipole shifting
    operator.

    i.e if we want the translated x-component of the dipole moment

    >>> from fmmgen.generator import generate_mappings
    >>> from fmmgen.expansions import M_shift
    >>> order = 2
    >>> x, y, z = sp.symbols('x y z')
    >>> idict, rdict = generate_mappings(order, [x, y, z])
    >>> M_shift((1, 0, 0), order, (x, y, z), idict)
    x*M[0, 0] + M[1, 0]
    """
    x, y, z = symbols
    modn = sum(n)

    grevlexkey = monomial_key('grevlex', symbols)
    monoms = sorted(itermonomials(symbols, modn), key=grevlexkey)
    result = sp.Integer(0)
    for monom in monoms:
        d = monom.as_powers_dict()
        k = d[x], d[y], d[z]
        nmink = n[0] - k[0], n[1] - k[1], n[2] - k[2]
        if nmink[0] >= 0 and nmink[1] >= 0 and nmink[2] >= 0:
            M = sp.MatrixSymbol('M', Nterms(order), 1)[index_dict[nmink]]
            sum_term = M * x**k[0] * y**k[1] * z**k[2] / \
                (factorial(k[0]) * factorial(k[1])*factorial(k[2]))
            result += sum_term
    return result


@functools.lru_cache(maxsize=None)
def Phi_derivatives(n, symbols):
    """
    Phi_derivatives(n, symbols)

    Returns symbolic expression for the nth derivative of 1/R

    Inputs:
    n, int:
        Tuple of derivatives
    symbols, tuple:
        Tuple of sympy symbols. Note that it *must* be a tuple for the
        functools.lru_cache to hash the input.

    Example:

    >>> import sympy as sp
    >>> from fmmgen.expansions import Phi_derivatives
    >>> x, y, z = sp.symbols('x y z')
    >>> Phi_derivatives((1, 0, 0), (x, y, z))
    -1.0*x/R**3
    """
    dx, dy, dz = symbols
    R = (dx**2 + dy**2 + dz**2)**(0.5)
    phi = 1/R
    deriv = sp.diff(phi, dx, n[0], dy, n[1], dz, n[2])
    deriv = deriv.subs(R, 'R')
    return deriv


def L(n, order, symbols, index_dict):
    dx, dy, dz = symbols
    modn = sum(n)
    modm_max = order - modn
    grevlexkey = monomial_key('grevlex', [dz, dy, dx])
    monoms = sorted(itermonomials(symbols, modm_max), key=grevlexkey)

    result = sp.Integer(0)
    for monom in monoms:
        d = monom.as_powers_dict()
        m = d[dx], d[dy], d[dz]
        npm = n[0] + m[0], n[1] + m[1], n[2] + m[2]
        if npm[0] >= 0 and npm[1] >= 0 and npm[2] >= 0:
            M = sp.MatrixSymbol('M', Nterms(order), 1)[index_dict[m]]
            result += M*Phi_derivatives(npm, symbols)

    return result


def L_shift(n, order, symbols, index_dict):
    x, y, z = symbols
    modn = sum(n)
    modk_max = order - modn
    grevlexkey = monomial_key('grevlex', [z, y, x])
    monoms = sorted(itermonomials([x, y, z], modk_max), key=grevlexkey)
    result = sp.Integer(0)
    for monom in monoms:
        d = monom.as_powers_dict()
        k = d[x], d[y], d[z]
        npk = n[0] + k[0], n[1] + k[1], n[2] + k[2]
        if npk[0] >= 0 and npk[1] >= 0 and npk[2] >= 0:
            L = sp.MatrixSymbol('L', Nterms(order), 1)[index_dict[npk]]
            sum_term = L * x**k[0] * y**k[1] * z**k[2] / \
                (factorial(k[0]) * factorial(k[1]) * factorial(k[2]))
            result += sum_term
    return result


def phi_deriv(order, symbols, index_dict, deriv=(0, 0, 0)):
    x, y, z = symbols
    grevlexkey = monomial_key('grevlex', [z, y, x])
    monoms = sorted(itermonomials([x, y, z], order), key=grevlexkey)
    monoms = itermonomials([x, y, z], order)

    result = sp.Integer(0)
    for monom in monoms:
        d = monom.as_powers_dict()
        n = d[x], d[y], d[z]
        L = sp.MatrixSymbol('L', Nterms(order), 1)[index_dict[n]]
        nmd = n[0] - deriv[0], n[1] - deriv[1], n[2] - deriv[2]
        if nmd[0] >= 0 and nmd[1] >= 0 and nmd[2] >= 0:
            sum_term = L * x**nmd[0] * y**nmd[1] * z**nmd[2] /  \
                      (factorial(nmd[0])*factorial(nmd[1])*factorial(nmd[2]))
            result += sum_term
    return result
