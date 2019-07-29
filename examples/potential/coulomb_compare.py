import numpy as np
import sys

n = int(sys.argv[1])

A = np.loadtxt(f'particles_n_{n}.txt', delimiter=',', usecols=(0, 1, 2, 3))


xv = A[:, 0]
yv = A[:, 1]
zv = A[:, 2]
qv = A[:, 3]

Phi = np.zeros(n)
for i in range(n):
    for j in range(n):
        if i != j:
            R = ((xv[i]-xv[j])**2 + (yv[i]-yv[j])**2 + (zv[i]-zv[j])**2)**0.5
            Phi[i] += qv[j]/R
print(Phi)
