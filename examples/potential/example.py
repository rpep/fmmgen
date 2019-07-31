import fmmgen

source_order = 1
order = source_order + 4
cse = False
atomic = True
precision='double'

fmmgen.generate_code(order, "operators",
                     precision=precision,
                     CSE=cse,
                     cython_wrapper=False,
                     potential=True,
                     field=True,
                     source_order=source_order,
                     atomic=atomic, minpow=5,
                     harmonic_derivs=True,
                     language='c++')
