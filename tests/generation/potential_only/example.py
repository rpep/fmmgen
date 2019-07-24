import fmmgen

order = 7

fmmgen.generate_code(order, "operators", CSE=False, generate_cython_wrapper=True)
