#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 08:41:10 2019

@author: ryan
"""
import os
import subprocess
import sympy as sp
# from utils import Nterms
from sympy.printing.ccode import C99CodePrinter
from sympy.printing.fcode import FCodePrinter
from sympy.printing.cxxcode import CXX17CodePrinter
from sympy.printing.llvmjitcode import LLVMJitPrinter
from utils import Nterms
import textwrap
from cse import cse
from generator import generate_mappings, generate_M_operators, \
                  generate_M_shift_operators, generate_L_operators, \
                  generate_L_shift_operators, \
                  generate_L2P_operators


q, x, y, z, R = sp.symbols('q x y z R')
symbols = (x, y, z)

language_mapping = {'c': C99CodePrinter,
                    'fortran': FCodePrinter,
                    'cxx': CXX17CodePrinter,
                    'llvm': LLVMJitPrinter
                    }


class FunctionPrinter:
    def __init__(self, language='c', precision='double', debug=True):
        self.debug = debug
        try:
            self.printer = language_mapping[language]()
        except KeyError:
            raise ValueError("Language not supported")

        self.precision = precision
        assert self.precision in ['float', 'double']

    def _generate_body(self, LHS, RHS, operator='='):
        # Note: replacement will *not* work with
        # common-subexpresion factoring.
        # Need to think more deeply about how
        # to do this. Given CSE seems to work less well
        # for expressions containing matrices, will
        # leave this til later to deal with.
        # sub_exprs, simp_rhs = sp.cse(RHS)
        #
        # code = self.printer.doprint(simp_RHS, assign_to=LHS)

        # Find the reduced RHS equation.
        code = ""

        if sp.symbols('R') in RHS.free_symbols:
            code += f'{self.precision} R = sqrt(x*x + y*y + z*z);\n'

        if not self.debug:
            sub_expressions, rRHS = cse(RHS)
            rRHS = sp.Matrix(rRHS)

            for var, _ in sub_expressions:
                code += f'{self.precision} {var};\n'

            for i, (var, sub_expr) in enumerate(sub_expressions):
                code += self.printer.doprint(sub_expr, assign_to=var) + "\n"
            code += self.printer.doprint(rRHS, assign_to=LHS).replace('=',
                                                                      operator)

        else:
            code += self.printer.doprint(RHS, assign_to=LHS).replace('=',
                                                                     operator)

        return code

    def _generate_header(self, name, LHS, RHS, inputs):
        types = []
        for arg in map(type, inputs):
            if arg == sp.MatrixSymbol:
                types.append(self.precision + ' *')
            else:
                types.append(self.precision)

        inputs.append(LHS)
        types.append(self.precision + ' *')

        combined_inputs = ', '.join([str(x) + ' ' + str(y) for x, y in
                                     zip(types, inputs)])
        return "void {}({})".format(name, combined_inputs)

    def generate(self, name, LHS, RHS, inputs, operator='='):
        header = self._generate_header(name, LHS, RHS, inputs)
        code = header + ' {\n'
        code += self._generate_body(LHS, RHS, operator)
        code += '\n}\n'

        header += ';\n'

        return header, code


def generate_code(order, name, generate_cython_wrapper=False, debug=True):
    p = FunctionPrinter(precision='float', debug=debug)

    header = ""
    body = ""

    x, y, z, q = sp.symbols('x y z q')
    symbols = (x, y, z)
    coords = [x, y, z]

    for i in range(order):
        idict, rdict = generate_mappings(i, symbols)

        M = sp.Matrix(generate_M_operators(i, symbols, idict))
        head, code = p.generate(f'P2M_{i}', 'M', M,
                                coords + [q], operator='+=')
        header += head
        body += code

        Ms = sp.Matrix(generate_M_shift_operators(i, symbols, idict))
        head, code = p.generate(f'M2M_{i}', 'Ms', Ms,
                                list(symbols) + \
                                [sp.MatrixSymbol('M', Nterms(i), 1)],
                                operator="+=")
        header += head
        body += code + '\n'

        L = sp.Matrix(generate_L_operators(i, symbols, idict))
        head, code = p.generate(f'M2L_{i}', 'L', L,
                               list(symbols) +  \
                               [sp.MatrixSymbol('M', Nterms(i), 1)],
                               operator="+=")
        header += head
        body += code + '\n'

        Ls = sp.Matrix(generate_L_shift_operators(i, symbols, idict))
        head, code = p.generate(f'L2L_{i}', 'Ls', Ls,
                               list(symbols) + \
                               [sp.MatrixSymbol('L', Nterms(i), 1)],
                               operator="+=")

        Fs = sp.Matrix(generate_L2P_operators(i, symbols, idict))
        head, code = p.generate(f'L2P_{i}', 'F', Fs,
                               list(symbols) + \
                               [sp.MatrixSymbol('L', Nterms(i), 1)],
                               operator="+=")

        header += head
        body += code + '\n'

    f = open(f"{name}.h", 'w')
    f.write("#ifndef FUNCTIONS_H\n#define FUNCTIONS_H\n")
    f.write(header)
    f.write("#endif")
    f.close()

    f = open(f"{name}.c", 'w')
    f.write("#include \"functions.h\"\n#include \"math.h\"\n")
    f.write(body)
    f.close()

    if generate_cython_wrapper:

        func_definitions = header.split(';\n')

        library = f"{name}"
        f = open(f"{name}wrap.pyx", "w")

        pyxcode = textwrap.dedent("""\
        # cython: language_level=3
        cimport numpy as np

        cdef extern from "{}.h":
            {}

        """).format(name, '\n    '.join(func_definitions))

        subsdict = {" *": "[:]",
                    "void": "cpdef",
                    "_": ""}

        for funcname in func_definitions:
            if not funcname:
                break
            pyfuncname = funcname
            for key, value in subsdict.items():
                pyfuncname = pyfuncname.replace(key, value)
            pyxcode += pyfuncname + ':\n'

            function_name = funcname.split("(")[0].split(" ")[1]
            args = funcname.split("(")[1][:-1].split(",")
            processed_args = []

            for arg in args:
                if "*" in arg:
                    arg = arg.replace("* ", "&") + "[0]"
                processed_args.append(arg.split(" ")[-1])

            pyxcode += '    ' + \
                       function_name + \
                       "(" + ', '.join(processed_args) + ')\n\n'

        f.write(pyxcode)
        f.close()


        print("writing build file")
        f = open(f"{name}wrap.pyxbld", "w")

        print(library)

        bldcode = textwrap.dedent("""\
        import numpy as np

        def make_ext(modname, pyxfilename):
            from distutils.extension import Extension
            return Extension(name = modname,
                             sources=[pyxfilename, '{}'],
                             include_dirs=[np.get_include(), '.'],
                             library_dirs=['.'],
                             extra_link_args=[],
                             extra_compile_args=['-O3'])
        """).format(library + '.c', library)

        f.write(bldcode)
        f.close()
