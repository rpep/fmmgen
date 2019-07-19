import numpy as np
import seaborn as sns
import sys
sns.set(style="darkgrid")
import matplotlib.pyplot as plt

column = int(sys.argv[1])

for i in range(1, 5+2, 2):
    V = np.loadtxt(f'error_order_{i}.txt', delimiter=',')
    sns.distplot(np.log10(np.sqrt(V[:, column]**2)), kde_kws=dict(cumulative=True), label=f'{i}')

plt.xlabel('$\log10(err)$')
plt.ylabel('Cumulative fraction of particles')
plt.legend(title='Expansion Order')
plt.show()
