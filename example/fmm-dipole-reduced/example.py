import fmmgen

order = 6
source_order = 1
cse = True
atomic = True
fmmgen.generate_code(order, "operators", CSE=cse,
                     generate_cython_wrapper=False,
                     potential=False,
                     field=True,
                     include_dir="include",
                     src_dir="src",
                     source_order=source_order,
                     atomic=atomic)
