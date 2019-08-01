import fmmgen.generator as gen
import fmmgen.expansions as exp
from fmmgen.utils import q, Nterms
import sympy as sp

x, y, z, R = sp.symbols('x y z R')
symbols = (x, y, z)

def test_L_shift_0_order_monopole_source():
    order = 0
    source = 0
    array_length = Nterms(order) - Nterms(source - 1)
    
    L = sp.MatrixSymbol('L', Nterms(order), 1)
    M_dict, _ = gen.generate_mappings(order, symbols, key='grevlex', source_order=source)
    L_dict, _ = gen.generate_mappings(order - source, symbols, key='grevlex', source_order=source)

    print(M_dict)
    print(L_dict)

    assert exp.L_shift((0, 0, 0), order, symbols, L_dict, source_order=source) == L[0]


def test_L_shift_1_order_monopole_source():
    order = 1
    source = 0
    array_length = Nterms(order) - Nterms(source - 1)

    L = sp.MatrixSymbol('L', Nterms(order), 1)
    M_dict, _ = gen.generate_mappings(order, symbols, key='grevlex', source_order=source)
    L_dict, _ = gen.generate_mappings(order - source, symbols, key='grevlex', source_order=source)

    print(M_dict)
    print(L_dict)

    assert exp.L_shift((0, 0, 0), order, symbols, L_dict, source_order=source) == L[0] + x*L[1] + y*L[2] + z*L[3]
    assert exp.L_shift((1, 0, 0), order, symbols, L_dict, source_order=source) == L[1]
    assert exp.L_shift((0, 1, 0), order, symbols, L_dict, source_order=source) == L[2]
    assert exp.L_shift((0, 0, 1), order, symbols, L_dict, source_order=source) == L[3]

def test_L_shift_1_order_dipole_source():
    order = 1
    source = 1
    array_length = Nterms(order) - Nterms(source - 1)

    L = sp.MatrixSymbol('L', Nterms(order), 1)
    M_dict, _ = gen.generate_mappings(order, symbols, key='grevlex', source_order=source)
    L_dict, _ = gen.generate_mappings(order - source, symbols, key='grevlex', source_order=0)

    print(f"M_dict = {M_dict}")
    print(f"L_dict = {L_dict}")

    assert exp.L_shift((0, 0, 0), order, symbols, L_dict, source_order=source) == L[0]
