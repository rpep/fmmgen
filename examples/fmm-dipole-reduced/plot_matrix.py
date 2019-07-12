import numpy as np
import matplotlib
import sys
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

n = int(sys.argv[1])
M = np.zeros((n, n))

m2l = np.loadtxt("m2lfile_10000.txt", delimiter=',', dtype=np.int)
p2p = np.loadtxt("p2pfile_10000.txt", delimiter=',', dtype=np.int)

for a, b in m2l:
    M[a, b] = 1

for a, b in p2p:
    M[a, b] = 2

plt.imshow(M)
plt.colorbar()
plt.show()
