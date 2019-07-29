import sympy as sp
from sympy.polys.orderings import monomial_key
from sympy import factorial
from sympy import itermonomials
from .utils import q, Nterms, generate_mappings
import functools


def M(n, symbols, source_order=0):
    """
    M(n, symbols)

    Returns symbolic expression for the n = (nx, ny, nz) multipole expansion.

    i.e if we want the x-component of the dipole moment:

    >>> x, y, z = sp.symbols('x y z')
    >>> M((1, 0, 0), (x, y, z))
    -q*x
    """
    dx, dy, dz = symbols

    if source_order > 0:
        raise NotImplementedError("Not implemented!")

    return q * (-1)**(n[0] + n[1] + n[2]) * dx**n[0] * dy**n[1] * dz**n[2]


def M_shift(n, order, symbols, index_dict, source_order=0):
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

    # Must iterate through the full set of monomial terms for the generation,
    # rather than using index_dict, because otherwise we miss terms.
    monoms, _ = generate_mappings(order, symbols, key='grevlex', source_order=0)

    x, y, z = symbols
    modn = sum(n)
    if modn < source_order:
        raise ValueError('sum(n) must be greater than or equal to source_order')
    result = sp.Integer(0)
    for i, k in enumerate(monoms.keys()):
        nmink = n[0] - k[0], n[1] - k[1], n[2] - k[2]
        if nmink[0] >= 0 and nmink[1] >= 0 and nmink[2] >= 0 and sum(nmink) >= source_order:
            array_index = index_dict[nmink]
            M = sp.MatrixSymbol('M', Nterms(order), 1)[array_index]
            sum_term = M * x**k[0] * y**k[1] * z**k[2] #/ \
#                (factorial(k[0]) * factorial(k[1])*factorial(k[2]))
            result += sum_term

    return result

def M_dipole(n, symbols, M_dict):
    x, y, z = symbols
    mux, muy, muz = sp.symbols('mux muy muz')
    order = max(sum(i) for i in M_dict)
    term = M_shift(n, order, symbols, M_dict, source_order=1)
    M = sp.MatrixSymbol('M', Nterms(order), 1)
    term = term.subs({M[M_dict[(1,0,0)]]: mux,
                         M[M_dict[(0,1,0)]]: muy,
                         M[M_dict[(0,0,1)]]: muz,
                        }
                       )

    # By default set to zero
    replacement_dict = {M[i]: 0 for i in range(Nterms(order))}
    # Replace terms with mux, muy, muz source terms.
    replacement_dict[M[M_dict[(1, 0, 0)]]] = mux
    replacement_dict[M[M_dict[(0, 1, 0)]]] = muy
    replacement_dict[M[M_dict[(0, 0, 1)]]] = muz
    # Coordinate transform to match P2M operators.
    replacement_dict[x] = -x
    replacement_dict[y] = -y
    replacement_dict[z] = -z
    return term.subs(replacement_dict)





@functools.lru_cache(maxsize=None)
def Phi_derivatives(n, symbols, harmonic=False):
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
    if not harmonic or n[2] < 2:
        dx, dy, dz = symbols
        R = (dx**2 + dy**2 + dz**2)**(0.5)
        phi = 1/R
        deriv = sp.diff(phi, dx, n[0], dy, n[1], dz, n[2])
        deriv = deriv.subs(R, 'R')
        return deriv
    else:
        k = (n[0], n[1], n[2] - 2)
        k1 = (k[0] + 2, k[1], k[2])
        k2 = (k[0], k[1] + 2, k[2])
        return -(Phi_derivatives(k1, symbols, harmonic=harmonic)) - (Phi_derivatives(k2, symbols, harmonic=harmonic))

def L(n, order, symbols, M_dict, eval_derivs=True, source_order=0, harmonic_derivs=False):
    assert order >= source_order, "order must be >= source_order"

    #print(sum(n), order - source_order)

    assert sum(n) <= order - source_order, "Terms are zero if sum(n) < order - source_order"

    dx, dy, dz = symbols
    modn = sum(n)
    modm_max = order - modn
    print(f'n = {n}')

    monoms, _ = generate_mappings(modm_max, symbols, key='grevlex', source_order=0)
    #print(monoms)

    if not eval_derivs:
        D = sp.MatrixSymbol('D', Nterms(order), 1)

    result = sp.Integer(0)
    for m in monoms.keys():
        npm = n[0] + m[0], n[1] + m[1], n[2] + m[2]
        #print(f'  m = {m}, npm = {npm}')
        if npm[0] >= 0 and npm[1] >= 0 and npm[2] >= 0 and sum(m) >= source_order:
            M = sp.MatrixSymbol('M', Nterms(order), 1)[M_dict[m]] / (factorial(m[0]) * factorial(m[1]) * factorial(m[2]))
            if eval_derivs:
                result += M*Phi_derivatives(npm, symbols, harmonic=harmonic_derivs)
            else:
                result += M*D[M_dict[npm]]
    return result


def L_shift(n, order, symbols, L_dict, source_order=0):
    x, y, z = symbols
    modn = sum(n)
    modk_max = order - modn

    # Must iterate through the full set of monomial terms for the generation,
    # rather than using index_dict, because otherwise we miss terms!
    monoms, _ = generate_mappings(modk_max, symbols, key='grevlex', source_order=0)

    # print(monoms.keys())

    # print(f"L_shift({n}), order = {order}, source_order={source_order}")
    result = sp.Integer(0)

    for k in monoms.keys():
        npk = n[0] + k[0], n[1] + k[1], n[2] + k[2]
        # print(f"  k = {k}, npk = {npk}")

        if npk[0] >= 0 and npk[1] >= 0 and npk[2] >= 0 and sum(npk) <= order-source_order:
            L = sp.MatrixSymbol('L', Nterms(order), 1)[L_dict[npk]]
            sum_term = L * x**k[0] * y**k[1] * z**k[2] / \
                (factorial(k[0]) * factorial(k[1]) * factorial(k[2]))
            result += sum_term
    return result


def phi_deriv(order, symbols, L_dict, deriv=(0, 0, 0), source_order=0):
    x, y, z = symbols
    # grevlexkey = monomial_key('grevlex', [z, y, x])
    # monoms = sorted(itermonomials([x, y, z], order), key=grevlexkey)
    # monoms = itermonomials([x, y, z], order)

    result = sp.Integer(0)
    for n in L_dict.keys():
        L = sp.MatrixSymbol('L', Nterms(order), 1)[L_dict[n]]
        nmd = n[0] - deriv[0], n[1] - deriv[1], n[2] - deriv[2]
        if nmd[0] >= 0 and nmd[1] >= 0 and nmd[2] >= 0:
            sum_term = L * x**nmd[0] * y**nmd[1] * z**nmd[2] /  \
                      (factorial(nmd[0])*factorial(nmd[1])*factorial(nmd[2]))
            result += sum_term
    return result
