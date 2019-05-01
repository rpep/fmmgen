import fmmgen.generator as gen
import fmmgen.expansions as exp
from fmmgen.utils import q, Nterms
import sympy as sp

x, y, z = sp.symbols('x y z')
symbols = (x, y, z)

def test_M_shift_monopole():
    order = 2
    source = 0
    array_length = Nterms(order) - Nterms(source - 1)
    
    M = sp.MatrixSymbol('M', Nterms(order), 1)
    M_dict, _ = gen.generate_mappings(order, symbols, key='grevlex', source_order=source)
    # Check monopole term
    assert exp.M_shift((0, 0, 0), order, symbols, M_dict, source_order=source) == M[0]
    # Check dipole terms
    assert exp.M_shift((1, 0, 0), order, symbols, M_dict, source_order=source) == M[1] + M[0]*x
    assert exp.M_shift((0, 1, 0), order, symbols, M_dict, source_order=source) == M[2] + M[0]*y
    assert exp.M_shift((0, 0, 1), order, symbols, M_dict, source_order=source) == M[3] + M[0]*z
    # Check quadrupole terms
    assert exp.M_shift((2, 0, 0), order, symbols, M_dict, source_order=source) ==  \
        x**2/2 * M[0] + x * M[1] + M[4]

    assert sp.expand(exp.M_shift((1, 1, 0), order, symbols, M_dict, source_order=source)) == x*y*M[0] +  y * M[1] + x * M[2] + M[5]

    assert sp.expand(exp.M_shift((1, 0, 1), order, symbols, M_dict, source_order=source)) == x*z*M[0] + z * M[1] + x*M[3] + M[6]
        
    

def test_M_shift_dipole():
    order = 2
    source = 1
    M = sp.MatrixSymbol('M', Nterms(order), 1)
    M_dict, _ = gen.generate_mappings(order, symbols, key='grevlex', source_order=source)

    # With a source order of 1, the monopole term is no longer produced. We should check
    # therefore, that M_shift returns an error.
    
    print(M_dict[(1, 0, 1)])

    try:
        exp.M_shift((0, 0, 0), order, symbols, M_dict, source_order=source)
        raise AssertionError("No exception raised!")
    except ValueError:
        pass
    
    assert exp.M_shift((1, 0, 0), order, symbols, M_dict, source_order=source) == M[0]
    assert exp.M_shift((0, 1, 0), order, symbols, M_dict, source_order=source) == M[1]
    assert exp.M_shift((0, 0, 1), order, symbols, M_dict, source_order=source) == M[2]

    assert sp.expand(exp.M_shift((2, 0, 0), order, symbols, M_dict, source_order=source)) == M[3] + M[0]*x

    assert sp.expand(exp.M_shift((1, 1, 0), order, symbols, M_dict, source_order=source)) == y * M[0] + x * M[1] + M[4]

    assert sp.expand(exp.M_shift((1, 0, 1), order, symbols, M_dict, source_order=source)) == z * M[0] + x*M[2] + M[5]
