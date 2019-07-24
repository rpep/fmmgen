import fmmgen

order = 7
cse = True

fmmgen.generate_code(order, "operators", CSE=cse, cython_wrapper=False,
                     include_dir="include", src_dir="src")
