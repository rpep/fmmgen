from  sympy.printing.ccode import C99CodePrinter as C99Base
from sympy.printing.cxxcode import CXX11CodePrinter as CXX11Base

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

language_mapping = {'c': CCodePrinter,
#                   'fortran': FCodePrinter,
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

    def _array(self, name, matrix, allocate=False, operator='=', atomic=False):
        code = ""

        if sp.symbols('R') in matrix.free_symbols:
            code += f'{self.precision} R = sqrt(x*x + y*y + z*z);\n'

        if allocate:
            code += f'{self.precision} {name}[{len(matrix)}];\n'

        if not self.debug:
            iterator = SymbolIterator(name)
            sub_expressions, rRHS = cse(matrix, symbols=iterator)
            rRHS = sp.Matrix(rRHS)
            for i, (var, sub_expr) in enumerate(sub_expressions):
                code += f'{self.precision} ' + self.printer.doprint(sub_expr, assign_to=var) + "\n"


            tmp = self.printer.doprint(rRHS, assign_to=name).replace('=',
                                                                    operator)

        else:
            tmp = self.printer.doprint(matrix, assign_to=name).replace('=',
                                                                       operator)

        if atomic:
            lines = tmp.split('\n')
            for l in lines:
                code += '#pragma omp atomic\n'
                code += l + '\n'
        else:
            code += tmp + '\n'
        return code




    def _generate_body(self, LHS, RHS, internal=[], operator='=', atomic=False):
        # Find the reduced RHS equation.
        logger.debug(f"Generating body for LHS = {str(LHS)}")
        code = ""

        for arr_name, matrix in internal:
            code += self._array(arr_name, matrix, allocate=True)

        code += self._array(LHS, RHS, operator=operator, atomic=atomic)
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

    def generate(self, name, LHS, RHS, inputs, operator='=', atomic=False, internal=[]):
        header = self._generate_header(name, LHS, RHS, inputs)
        code = header + ' {\n'
        code += self._generate_body(LHS, RHS, internal, operator, atomic=atomic)
        code += '\n}\n'
        header += ';\n'

        return header, code
