import numpy as np
import matplotlib.pyplot as plt
f = open("output.txt")
a = f.read().split('\n')

Nparticles = []
t_direct = []
t_approx = []

for l in a:
    if 'Nparticles' in l:
        Nparticles.append(int(l.split(' ')[-1]))
    elif 't_direct' in l:
        t_direct.append(float(l.split(' ')[-1]))
    elif 'Approx.' in l:
        print(l)
        t_approx.append(float(l.split(' ')[-3]))



plt.loglog(Nparticles, t_direct, marker='o', label='direct')
plt.loglog(Nparticles, t_approx, marker='o', label='approx')
plt.show()
