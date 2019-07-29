#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  9 08:41:10 2019

@author: ryan
"""
import os
import subprocess
import sympy as sp
import logging
# from sympy.printing.fcode import FCodePrinter
from fmmgen.printers import CCodePrinter, CXXCodePrinter
from fmmgen.utils import Nterms
import textwrap
from fmmgen.cse import cse
from fmmgen.generator import generate_mappings, generate_M_operators, \
                  generate_M_shift_operators, generate_L_operators, \
                  generate_L_shift_operators, \
                  generate_L2P_operators, \
                  generate_M2P_operators, \
                  generate_P2P_operators


logger = logging.getLogger(name="fmmgen")


q, x, y, z, R = sp.symbols('q x y z R')
symbols = (x, y, z)

language_mapping = {'c': CCodePrinter,
#                   'fortran': FCodePrinter,
                    'c++': CXXCodePrinter,
                    }



class FunctionPrinter:
    def __init__(self, language='c', precision='double', debug=True, gpu=False, minpow=False):
        logger.info(f"Function Printer created with precision \"{precision}\"")

        self.gpu = gpu
        if self.gpu:
            logger.info(f"Writing CUDA __device__ functions is enabled")

        if not debug:
            logger.info(f"CSE is enabled")
        else:
            logger.info(f"CSE is disabled")
        self.debug = debug

        try:
            if minpow:
                self.printer = language_mapping[language](minpow=minpow)
            else:
                self.printer = language_mapping[language]()
        except KeyError:
            raise ValueError("Language not supported")

        self.precision = precision
        assert self.precision in ['float', 'double']

    def _generate_body(self, LHS, RHS, operator='=', atomic=False):
        # Find the reduced RHS equation.
        logger.debug(f"Generating body for LHS = {str(LHS)}")
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


            tmp = self.printer.doprint(rRHS, assign_to=LHS).replace('=',
                                                                      operator)
        else:
            tmp = self.printer.doprint(RHS, assign_to=LHS).replace('=',
                                                                     operator)
        if atomic:
            lines = tmp.split('\n')
            for l in lines:
                code += '#pragma omp atomic\n'
                code += l + '\n'
        else:
            code += tmp

        return code

    def _generate_header(self, name, LHS, RHS, inputs):
        logger.debug(f"Generating headerfile for LHS = {str(LHS)}")
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

        if self.gpu:
            return "__device__ void {}({})".format(name, combined_inputs)
        else:
            return "void {}({})".format(name, combined_inputs)

    def generate(self, name, LHS, RHS, inputs, operator='=', atomic=False):
        header = self._generate_header(name, LHS, RHS, inputs)
        code = header + ' {\n'
        code += self._generate_body(LHS, RHS, operator, atomic=atomic)
        code += '\n}\n'
        header += ';\n'

        return header, code


def generate_code(order, name, precision='double',
                  cython_wrapper=False,
                  CSE=False, harmonic_derivs=False,
                  include_dir=None, src_dir=None,
                  potential=True, field=True,
                  source_order=0, atomic=False,
                  gpu=False, minpow=0, language='c'):
    """
    Inputs:

    order, int:
        Expansion order of the FMM operators.

    name, str:
        Name of generated files.

        Note: Cython compiles *.pyx files to the *.c files with the same name,
        which means that we have to use a different name for the cython file.
        We therefore append "wrap" to the end of the name for the cython file,
        and this is therefore the name of the Python module which must be
        imported if using pyximport.

    cython_wrapper, bool:
        Enable generation of a Cython wrapper for the C files.

    CSE, bool:
        Enable common subexpression elimination, to reduce the op count in
        the generated code.

    harmonic_derivs, bool:
        The harmonicity of the Laplace means that at order p, there are only
        2p - 1 independent derivatives. Enabling this option therefore computes
        some derivatives as combinations of others.

    source_order, int:
        If source_order > 0 then we set certain multipole terms to zero
        in the local expansion, and hence they are not used. This is useful if,
        for example, we only have pure dipoles or quadrupoles in the system.

    minpow, int:
        If minpow is set, pow(x, n) expressions are expanded such that
        if n < minpow, the expression in code is printed as multiplications.

        e.g. if a sympy expression is pow(x, 2) + pow(y, 6) and minpow is 5,
        the printed version will be x*x + pow(y, 6)

    """
    assert language in ['c', 'c++'], "Language must be 'c' or 'c++'"
    if language == 'c':
        fext = 'c'
        hext = 'h'
    if language == 'c++':
        fext = 'cpp'
        hext = 'h'
    
    
    logger.info(f"Generating FMM operators to order {order}")
    assert precision in ['double', 'float'], "Precision must be float or double"
    logger.info(f"Precision = {precision}")
    if CSE:
        logger.info(f"CSE Enabled")
        p = FunctionPrinter(precision=precision, debug=False, minpow=minpow)
    else:
        logger.info(f"CSE Disabled")
        p = FunctionPrinter(precision=precision, debug=True, minpow=minpow)

    header = ""
    body = ""

    x, y, z, q = sp.symbols('x y z q')
    symbols = (x, y, z)
    coords = [x, y, z]


    start = source_order
    if field and not potential:
        # No point starting at source_order
        # because no field calculation can be done
        # at this multipole order - the L2P derivative
        # is 0.
        start += 1

    for i in range(start, order):
        logger.info(f"Generating order {i} operators")
        M_dict, _ = generate_mappings(i, symbols,'grevlex',  source_order=source_order)
        L_dict, _ = generate_mappings(i - source_order, symbols, 'grevlex', source_order=0)



        M = sp.Matrix(generate_M_operators(i, symbols, M_dict))

        # print(f"M = {M}")

        head, code = p.generate(f'P2M_{i}', 'M', M,
                                coords + [q], operator='+=')
        header += head
        body += code

        Ms = sp.Matrix(generate_M_shift_operators(i, symbols, M_dict, source_order=source_order))

        # print(f"Ms = {Ms}")

        head, code = p.generate(f'M2M_{i}', 'Ms', Ms,
                                list(symbols) + \
                                [sp.MatrixSymbol('M', Nterms(i), 1)],
                                operator="+=", atomic=atomic)
        header += head
        body += code + '\n'

        L = sp.Matrix(generate_L_operators(i, symbols, M_dict, L_dict,
                      source_order=source_order, harmonic_derivs=harmonic_derivs))

        # print(f"L = {L}")


        head, code = p.generate(f'M2L_{i}', 'L', L,
                               list(symbols) +  \
                               [sp.MatrixSymbol('M', Nterms(i), 1)],
                                operator="+=", atomic=atomic)
        header += head
        body += code + '\n'


        Ls = sp.Matrix(generate_L_shift_operators(i, symbols, L_dict, source_order=source_order))
        head, code = p.generate(f'L2L_{i}', 'Ls', Ls,
                               list(symbols) + \
                               [sp.MatrixSymbol('L', Nterms(i), 1)],
                                operator="+=", atomic=atomic)

        header += head
        body += code + '\n'

        L2P = generate_L2P_operators(i, symbols, L_dict,
                                     potential=potential,
                                     field=field)

        Fs = sp.Matrix(L2P)
        head, code = p.generate(f'L2P_{i}', 'F', Fs,
                               list(symbols) + \
                               [sp.MatrixSymbol('L', Nterms(i), 1)],
                                operator="+=", atomic=atomic)

        header += head
        body += code + '\n'

        M2P = generate_M2P_operators(i, symbols, M_dict,
                                     potential=potential,
                                     field=field, source_order=source_order, harmonic_derivs=harmonic_derivs)
        Fs = sp.Matrix(M2P)
        head, code = p.generate(f'M2P_{i}', 'F', Fs,
                                list(symbols) + \
                                [sp.MatrixSymbol('M', Nterms(i), 1)],
                                operator="+=", atomic=atomic)
        header += head
        body += code + '\n'

        if i == start:
            P2P = sp.Matrix(generate_P2P_operators(symbols, M_dict,
                                                   potential=potential,
                                                   field=field,
                                                   source_order=source_order))

            head, code = p.generate(f'P2P', 'F', P2P,
                                    list(symbols) + \
                                    [sp.MatrixSymbol('S', Nterms(i), 1)],
                                    operator="+=", atomic=atomic
                                    )
            header += head
            body += code + '\n'

    # print("Here!")
    # We now generate wrapper functions that cover all orders generated.
    unique_funcs = []
    func_definitions = header.split(';\n')
    # print(f"func_defs = {func_definitions}")
    for func in func_definitions:
        # Must do it this way in order to avoid breaking
        # for expansions > 10.
        function_name = func.split('(')[0]
        print(f"Function_name = {function_name}")
        end_string = f'_{start}'
        if end_string == function_name[-len(end_string):]:
            print("Unique!")
            unique_funcs.append(func)
        else:
            print(f"{func} not unique")
            print(f"  {end_string}  {function_name[-len(end_string)-1:]}")

    wrapper_funcs = [f.replace(')', ', int order)').replace(f'_{start}', '')
                     for f in unique_funcs]

    print(wrapper_funcs)

    func_definitions += wrapper_funcs
    print('\n'.join(func_definitions))

    for wfunc, func in zip(wrapper_funcs, unique_funcs):
        # Add to header file
        header += wfunc + ';\n'
        # Create a switch statement that covers all functions:
        code = wfunc + " {\n"
        code += 'switch (order) {\n'
        for i in range(start, order):
            code += '  case {}:\n'.format(i)
            print(func)
            replaced_code = func.replace(f'_{start}', f'_{i}').replace('* ','').replace('double ','').replace('float ','').replace('void ', '')
            print(f"replaced_code: {replaced_code}")
            code += '    ' + replaced_code + ';\n    break;\n'
        code += "  }\n}\n"
        # print(code)
        body += code

    if not include_dir:
        f = open(f"{name}.{hext}", 'w')
    else:
        f = open(f"{include_dir.rstrip('/')}/{name}.{hext}", 'w')
    f.write(f"#pragma once\n")
    f.write(f"#define FMMGEN_MINORDER {start}\n")
    f.write(f"#define FMMGEN_MAXORDER {order}\n")
    f.write(f"#define FMMGEN_SOURCE_SIZE {Nterms(source_order) - Nterms(source_order - 1)}\n")
    if potential and not field:
        osize = 1
    elif field and not potential:
        osize = 3
    elif field and potential:
        osize = 4
    f.write(f"#define FMMGEN_OUTPUT_SIZE {osize}\n")
    f.write(header)
    f.close()

    if not src_dir:
        f = open(f"{name}.{fext}", 'w')
    else:
        f = open(f"{src_dir.rstrip('/')}/{name}.{fext}", 'w')

    f.write(f'#include "{name}.{hext}"\n')
    if language == 'c':
        f.write(f'#include "math.h"\n')
    elif language == 'c++':
        f.write(f'#include<cmath>\n')
    
    f.write(body)
    f.close()

    if cython_wrapper and gpu:
        raise Warning("Cannot write a Cython wrapper for GPU code; skipping")

    elif cython_wrapper:
        logger.info(f"Generating Cython wrapper: {name}_wrap.pyx")
        library = f"{name}"

        f = open(f"{name}_decl.pxd", "w")
        pxdcode = textwrap.dedent("""\
        cdef extern from "{}.h":
            cdef int FMMGEN_MINORDER
            cdef int FMMGEN_MAXORDER
            cdef int FMMGEN_SOURCE_SIZE
            cdef int FMMGEN_OUTPUT_SIZE
            {}
        """)
        f.write(pxdcode.format(name, '\n    '.join(func_definitions)))
        
        f.close()

        f = open(f"{name}_wrap.pyx", "w")
        # expose the C functions from the header file.
        pyxcode = textwrap.dedent("""\
        # cython: language_level=3
        cimport numpy as np
        cimport {}

        FMMGEN_MINORDER = {}.FMMGEN_MINORDER
        FMMGEN_MAXORDER = {}.FMMGEN_MAXORDER
        FMMGEN_SOURCE_SIZE = {}.FMMGEN_SOURCE_SIZE
        FMMGEN_OUTPUT_SIZE = {}.FMMGEN_OUTPUT_SIZE
        """).format(*[name + '_decl']*5) 

        subsdict = {" *": "[:]",
                    "void": "cpdef",
                    "_": ""}

        # Generate the actual wrapper code
        for funcname in func_definitions:
            # print(funcname)
            if not funcname:
                continue
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
                       name + '_decl.' + function_name + \
                       "(" + ', '.join(processed_args) + ')\n\n'

        f.write(pyxcode)
        f.close()

        f = open(f"{name}_wrap.pyxbld", "w")

        # print(library)

        logger.info(f"Generating Cython buildfile: {name}_wrap.pyxbld")
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
