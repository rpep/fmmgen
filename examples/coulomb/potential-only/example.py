import fmmgen

order = 6
source_order = 0
cse = True
atomic = True
precision='double'

fmmgen.generate_code(order, "operators",
                     precision=precision,
                     CSE=cse,
                     cython_wrapper=False,
                     potential=True,
                     field=False,
                     source_order=source_order,
                     atomic=atomic, minpow=5,
                     harmonic_derivs=True,
                     language='c++')
