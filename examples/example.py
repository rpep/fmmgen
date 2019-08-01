import fmmgen

source_order = 2
order = source_order + 8
cse = True
atomic = True
precision='double'

fmmgen.generate_code(order, "operators",
                     precision=precision,
                     CSE=cse,
                     cython=False,
                     potential=True,
                     field=True,
                     source_order=source_order,
                     atomic=atomic, minpow=5,
                     harmonic_derivs=True,
                     language='c++')
