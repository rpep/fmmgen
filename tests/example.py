import fmmgen

order = 3

fmmgen.generate_code(order, "operators", CSE=True, generate_cython_wrapper=True)
