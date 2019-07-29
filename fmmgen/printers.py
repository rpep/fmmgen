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
