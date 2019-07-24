import numpy as np
import seaborn as sns
import sys
sns.set(style="darkgrid")
import matplotlib.pyplot as plt


print(sys.argv)

n = int(sys.argv[1])
ncrit = int(sys.argv[2])
theta = float(sys.argv[3])
maxorder = int(sys.argv[4])
column = int(sys.argv[5])

for i in range(2, maxorder, 2):
    try:
        filename = f'errors_lazy_p_{i}_n_{n}_ncrit_{ncrit}_theta_{theta:.6f}.txt'
        print(filename)
        V = np.loadtxt(filename, delimiter=',')
        sns.distplot(np.log10(np.sqrt(V[:, column]**2)), kde_kws=dict(cumulative=True))
        print(f'Plotted {i}')
    except:
        print(f'Failed {i}')
        pass

#plt.xlim("$log(\epsilon)$")
#plt.ylim("Cumulative fraction of particles")
plt.show()
