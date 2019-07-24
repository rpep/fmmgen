import fmmgen

order = 7
source_order = 1
cse = True
atomic = True
precision='double'

fmmgen.generate_code(order, "operators",
                     precision=precision,
                     CSE=cse,
                     cython_wrapper=False,
                     potential=False,
                     field=True,
                     include_dir="include",
                     src_dir="src",
                     source_order=source_order,
                     atomic=atomic, minpow=4)
