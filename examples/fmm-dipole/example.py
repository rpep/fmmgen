import fmmgen

order = 8
cse = False

fmmgen.generate_code(order, "operators", CSE=cse, cython_wrapper=False,
                     include_dir="include", src_dir="src")
