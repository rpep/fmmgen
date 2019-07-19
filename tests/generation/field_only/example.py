import fmmgen

order = 5

fmmgen.generate_code(order, "operators", CSE=False, generate_cython_wrapper=True, potential=False, field=True)

