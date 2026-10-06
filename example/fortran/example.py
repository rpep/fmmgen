import fmmgen

source_order = 0
order = source_order + 12
cse = True
precision = 'double'

# compress and planar are both needed: the driver picks the variant at run time
# (--compress, --dim), exactly as the C++ example does.
fmmgen.generate_code(order, "operators",
                     precision=precision,
                     CSE=cse,
                     potential=True,
                     field=True,
                     source_order=source_order,
                     minpow=11,
                     harmonic_derivs=True,
                     compress=True,
                     planar=True,
                     language='fortran')
