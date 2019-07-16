from  sympy.printing.ccode import C99CodePrinter as C99Base

class C99CodePrinter(C99Base):
    def __init__(self, settings={}, minpow=False):
        super(C99CodePrinter, self).__init__(settings)
        self.minpow = minpow
        
    def _print_Pow(self, expr):
        if self.minpow:
            if expr.exp.is_integer and expr.exp > 0 and expr.exp <= self.minpow:
                return '(' + '*'.join([self._print(expr.base) for i in range(expr.exp)]) + ')'
            else:
                return super()._print_Pow(expr)
        else:
            return super()._print_Pow(expr)
