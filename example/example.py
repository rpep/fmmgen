import fmmgen

source_order = 0
order = source_order + 5
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
                     harmonic_derivs=False,
                     language='c++')
