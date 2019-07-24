import fmmgen

order = 5
source_order = 1
cse = True
atomic = True
precision='double'
potential=False
field=True
cython_wrapper=False

fmmgen.generate_code(order, "operators",
                     precision=precision,
                     CSE=cse,
                     cython_wrapper=cython_wrapper,
                     potential=potential,
                     field=field,
                     include_dir="include",
                     src_dir="src",
                     source_order=source_order,
                     atomic=atomic)
