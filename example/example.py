import fmmgen

source_order = 0
order = source_order + 12
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
                     atomic=atomic, minpow=11,
                     harmonic_derivs=True,
                     language='c++')
