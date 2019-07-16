import fmmgen

order = 4
cse = True
atomic = False

fmmgen.generate_code(order, "operators", CSE=cse, generate_cython_wrapper=False,
                     include_dir="include", src_dir="src", atomic=atomic)
