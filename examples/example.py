import os
import fmmgen

source_order = 0
order = source_order + 12
cse = True
atomic = False
precision='double'

fmmgen.generate_code(order, "operators",
                     precision=precision,
                     CSE=cse,
                     cython=False,
                     potential=True,
                     field=True,
                     source_order=source_order,
                     atomic=atomic, minpow=11,
                     harmonic_derivs=True,
                     compress=True,
                     planar=True,
                     language='c++',
                     # Speed knobs (output is identical): jobs = worker processes (0 = all
                     # CPUs); split = one operators_<operator>.cpp file per operator so `make -j`
                     # can compile them in parallel (False = one operators.cpp, as before).
                     jobs=int(os.environ.get('FMMGEN_JOBS', 1)),
                     split=os.environ.get('FMMGEN_SPLIT', '0') == '1',
                     cache_dir=os.environ.get('FMMGEN_CACHE', 'auto') or None)
