from sympy import itermonomials
from .expansions import M, M_shift, L, L_shift, phi_deriv
from sympy.polys.orderings import monomial_key

def generate_mappings(order, symbols, key='grevlex'):
    x, y, z = symbols
    monoms = itermonomials([x, y, z], order)
    if key:
        monom_key = monomial_key(key, [z, y, x])
        monoms = sorted(monoms, key=monom_key)

    index_dict = {}
    rindex_dict = {}
    for i, monom in enumerate(monoms):
        d = monom.as_powers_dict()
        n = d[x], d[y], d[z]
        index_dict[n] = i
        rindex_dict[i] = n
    return index_dict, rindex_dict


def generate_M_operators(order, symbols, index_dict):
    x, y, z = symbols
    M_operators = []
    for n in index_dict.keys():
        M_operators.append(M(n, symbols))
    return M_operators


def generate_M_shift_operators(order, symbols, index_dict):
    x, y, z = symbols
    M_operators = []
    for n in index_dict.keys():
        M_operators.append(M_shift(n, order, symbols, index_dict))
    return M_operators


def generate_L_operators(order, symbols, index_dict):
    x, y, z = symbols
    L_operators = []
    for n in index_dict.keys():
        L_operators.append(L(n, order, symbols, index_dict))
    return L_operators


def generate_L_shift_operators(order, symbols, index_dict):
    x, y, z = symbols
    L_shift_operators = []
    for n in index_dict.keys():
        L_shift_operators.append(L_shift(n, order, symbols, index_dict))
    return L_shift_operators


def generate_L2P_operators(order, symbols, index_dict):
    phi = phi_deriv(order, symbols, index_dict, deriv=(0, 0, 0))
    Fx = -phi_deriv(order, symbols, index_dict, deriv=(1, 0, 0))
    Fy = -phi_deriv(order, symbols, index_dict, deriv=(0, 1, 0))
    Fz = -phi_deriv(order, symbols, index_dict, deriv=(0, 0, 1))
    return [phi, Fx, Fy, Fz]
