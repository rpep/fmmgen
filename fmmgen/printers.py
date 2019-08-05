
from sympy.printing.ccode import C99CodePrinter as C99Base
from sympy.printing.cxxcode import CXX11CodePrinter as CXX11Base
import logging
logger = logging.getLogger(name="fmmgen")
import sympy as sp
from fmmgen.cse import cse
from fmmgen.opts import basic as opts
from sympy import count_ops
from sympy.printing.precedence import precedence
from sympy.core.mul import _keep_coeff

class CCodePrinter(C99Base):
    def __init__(self, settings={}, minpow=False):
        super(C99Base, self).__init__(settings)
        self.minpow = minpow

    def _print_Pow(self, expr):
        if self.minpow:
            if expr.exp.is_integer and expr.exp > 0 and expr.exp <= self.minpow:
                return '(' + '*'.join([self._print(expr.base) for i in range(expr.exp)]) + ')'

            elif expr.exp.is_integer and expr.exp < 0 and expr.exp >= -self.minpow:
                expr = '(1 / (' + '*'.join([self._print(expr.base) for i in range(abs(expr.exp))]) + '))'
                return expr
            else:
                return super()._print_Pow(expr)
        else:
            return super()._print_Pow(expr)

    def _print_Mul(self, expr):
        prec = precedence(expr)
        c, e = expr.as_coeff_Mul()

        if c == 1.0:
            expr = e
            sign = ""
        elif e == 1.0:
            expr = c
            sign = ""
        
        elif c < 0:
            if c == -1.0:
                expr = e
                sign = "-"
            elif e == -1.0:
                expr = c
                sign = "-"
            else:
                expr = _keep_coeff(-c, e)
                sign = "-"
        else:
            sign = ""

        a = []  # items in the numerator
        b = []  # items that are in the denominator (if any)

        pow_paren = []  # Will collect all pow with more than one base element and exp = -1

        if self.order not in ('old', 'none'):
            args = expr.as_ordered_factors()
        else:
            # use make_args in case expr was something like -x -> x
            args = Mul.make_args(expr)

        # Gather args for numerator/denominator
        for item in args:
            if item.is_commutative and item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
                if item.exp != -1:
                    b.append(Pow(item.base, -item.exp, evaluate=False))
                else:
                    if len(item.args[0].args) != 1 and isinstance(item.base, Mul):   # To avoid situations like #14160
                        pow_paren.append(item)
                    b.append(Pow(item.base, -item.exp))
            else:
                a.append(item)

        a = a or [S.One]

        a_str = [self.parenthesize(x, prec) for x in a]
        b_str = [self.parenthesize(x, prec) for x in b]

        # To parenthesize Pow with exp = -1 and having more than one Symbol
        for item in pow_paren:
            if item.base in b:
                b_str[b.index(item.base)] = "(%s)" % b_str[b.index(item.base)]

        if not b:
            return sign + '*'.join(a_str)
        elif len(b) == 1:
            return sign + '*'.join(a_str) + "/" + b_str[0]
        else:
            return sign + '*'.join(a_str) + "/(%s)" % '*'.join(b_str)


class CXXCodePrinter(CXX11Base):
    def __init__(self, settings={}, minpow=False):
        super(CXX11Base, self).__init__(settings)
        self.minpow = minpow

    def _print_Pow(self, expr):
        if self.minpow:
            if expr.exp.is_integer and expr.exp > 0 and expr.exp <= self.minpow:
                return '(' + '*'.join([self._print(expr.base) for i in range(expr.exp)]) + ')'

            elif expr.exp.is_integer and expr.exp < 0 and expr.exp >= -self.minpow:
                expr = '(1 / (' + '*'.join([self._print(expr.base) for i in range(abs(expr.exp))]) + '))'
                return expr
            else:
                return super()._print_Pow(expr)
        else:
            return super()._print_Pow(expr)

    def _print_Mul(self, expr):
        prec = precedence(expr)
        c, e = expr.as_coeff_Mul()

        if c == 1.0:
            return str(e)
        elif e == 1.0:
            return str(c)
        
        if c < 0:
            expr = _keep_coeff(-c, e)
            sign = "-"
        else:
            sign = ""

        a = []  # items in the numerator
        b = []  # items that are in the denominator (if any)

        pow_paren = []  # Will collect all pow with more than one base element and exp = -1

        if self.order not in ('old', 'none'):
            args = expr.as_ordered_factors()
        else:
            # use make_args in case expr was something like -x -> x
            args = Mul.make_args(expr)

        # Gather args for numerator/denominator
        for item in args:
            if item.is_commutative and item.is_Pow and item.exp.is_Rational and item.exp.is_negative:
                if item.exp != -1:
                    b.append(Pow(item.base, -item.exp, evaluate=False))
                else:
                    if len(item.args[0].args) != 1 and isinstance(item.base, Mul):   # To avoid situations like #14160
                        pow_paren.append(item)
                    b.append(Pow(item.base, -item.exp))
            else:
                a.append(item)

        a = a or [S.One]

        a_str = [self.parenthesize(x, prec) for x in a]
        b_str = [self.parenthesize(x, prec) for x in b]

        # To parenthesize Pow with exp = -1 and having more than one Symbol
        for item in pow_paren:
            if item.base in b:
                b_str[b.index(item.base)] = "(%s)" % b_str[b.index(item.base)]

        if not b:
            return sign + '*'.join(a_str)
        elif len(b) == 1:
            return sign + '*'.join(a_str) + "/" + b_str[0]
        else:
            return sign + '*'.join(a_str) + "/(%s)" % '*'.join(b_str)

language_mapping = {'c': CCodePrinter,
                    'c++': CXXCodePrinter,
                    }



class SymbolIterator:
    def __init__(self, name):
        self.name = name
        self.num = 0

    def __iter__(self):
        return self

    def __next__(self):
        num = self.num
        self.num += 1
        return sp.Symbol(self.name + 'tmp' + str(num))


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

    def _array(self, name, matrix, allocate=False, operator='=', atomic=False, ignore_symbols=[]):
        opscount = 0
        code = ""

        
        # Testing on Godbolt with GCC 9.1 and ICPC shows that
        # pow(x, 0.5) generates fewer instructions than
        # sqrt(x), so will swap. However, leave R here in case
        # it gets used in expansions in future for some reason.

        light_ignore=[]
        if sp.symbols('R') in matrix.free_symbols:
            code += f'{self.precision} R = sqrt(x*x + y*y + z*z);\n'
            light_ignore.append('R')
        if sp.symbols('Rinv') in matrix.free_symbols:
            code += f'{self.precision} Rinv = pow(x*x + y*y + z*z, -0.5);\n'
            light_ignore.append('Rinv')
            
        if allocate:
            code += f'{self.precision} {name}[{len(matrix)}];\n'

        if not self.debug:
            # print('Printing with CSE')
            iterator = SymbolIterator(name)
            # print(f'ignoring {name} in cse')
            sub_expressions, rmatrix = cse(matrix, optimizations=opts,
                                           symbols=iterator,
                                           ignore=(ignore_symbols),
                                           light_ignore=light_ignore)
        
            rmatrix = sp.Matrix(rmatrix)
            for i, (var, sub_expr) in enumerate(sub_expressions):
                opscount += count_ops(sub_expr)
                code += f'{self.precision} ' + self.printer.doprint(sub_expr, assign_to=var) + "\n"

            opscount += count_ops(rmatrix)
            tmp = self.printer.doprint(rmatrix, assign_to=name).replace('=',
                                                                        operator)

        else:
            # print('Printing without CSE')
            opscount += count_ops(matrix)
            tmp = self.printer.doprint(matrix, assign_to=name).replace('=',
                                                                       operator)

        if atomic:
            lines = tmp.split('\n')
            for l in lines:
                code += '#pragma omp atomic\n'
                code += l + '\n'
        else:
            code += tmp + '\n'
        return code, opscount




    def _generate_body(self, LHS, RHS, internal=[], operator='=', atomic=False, ignore=[]):
        # Find the reduced RHS equation.
        opscount = 0
        logger.debug(f"Generating body for LHS = {str(LHS)}")
        code = ""
        
        for arr_name, matrix in internal:
            codetext, ops = self._array(arr_name, matrix, allocate=True, ignore_symbols=[arr_name]+ignore)
            code += codetext
            opscount += ops
            
        codetext, ops = self._array(LHS, RHS, operator=operator, atomic=atomic)
        code += codetext
        opscount += ops
        return code, opscount


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

    def generate(self, name, LHS, RHS, inputs, operator='=', atomic=False, internal=[], ignore=[]):
        header = self._generate_header(name, LHS, RHS, inputs)
        code = header + ' {\n'
        codetext, opscount = self._generate_body(LHS, RHS, internal, operator, atomic=atomic, ignore=ignore)
        code += codetext
        code += '\n}\n'
        header += ';\n'

        return header, code, opscount
