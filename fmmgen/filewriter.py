#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 15:50:03 2019

@author: ryan
"""
import sympy as sp
from utils import Nterms
from generator import generate_mappings, generate_M_operators, \
                      generate_M_shift_operators, generate_L_operators, \
                      generate_L_shift_operators, generate_L2P_operators


extensions = {"c": ("c", "h"),
              "cxx": ("cpp", "hpp"),
              "fortran": ("f90", "h90")}


class FileWriter:
    def __init__(self, symbols, order=2, language='c'):
        self.order = order
        self.index_dict, self.rindex_dict = generate_mappings(order, symbols)
        self.printer = FunctionPrinter(language=language)
        self.f_ext, self.h_ext = extensions[language]

    def _generate_code_func(self, func, out_label, func_basename):
        n = Nterms(order)
        header = ""
        code = ""
        for i in range(self.order):
            arr = sp.Matrix(func(order, symbols, index_dict))
            h, c = self.printer.generate(f'{func_basename}_{order}', out_label, arr)
            header += h
            code += c

        return header, code


    def generate_code(self):
        header = ""
        code = ""

        funcs = [generate_M_operators, generate_M_shift_operators,
                 generate_L_operators, generate_L_shift_operators,
                 generate_L2P_operators]

        labels = 'M', 'M', 'L', 'L', 'F'
        basenames = 'P2M', 'M2M', 'M2L', 'L2L', 'L2P'

        for func, lab, base in zip(funcs, labels, basenames):
            h, c = self._generate_code_func(func, lab, base)
            header += h
            code += c

        print(header)
        print(code)
