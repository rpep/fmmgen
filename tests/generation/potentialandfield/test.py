import pyximport
pyximport.install()
import numpy as np
import operators_wrap as fmm
from fmmgen.utils import Nterms



order = 5
x = 1.0
y = 2.0
z = 5.0
q = 3.0
M = np.zeros(Nterms(order))

fmm.P2M(x, y, z, q, M, order)

print(M)
