from sympy.simplify import powsimp
from sympy.simplify.cse_opts import sub_pre, sub_post
from sympy.core.exprtools import factor_terms
from sympy.polys.polyfuncs import horner

# By default, CSE can miss some obvious optimisations
# due to how complex some of the expressions can be.
# By using preprocessing, we can simplify the expressions
# in some ways. We don't want to use sp.simplify, because
# these optimisations are too broad (i.e. we have no trig
# functions, so there's no point passing over with these)
# ; we can use instead a select few.

# powsimp changes things like (x^2)^4 into x^8
# ratsimp
basic = [
    # (horner, None),
    (powsimp, None),
    (sub_pre, sub_post),
    (factor_terms, None),
]
