import numpy as np
import seaborn as sns
import sys
sns.set(style="darkgrid")
import matplotlib.pyplot as plt

column = int(sys.argv[1])

for i in range(2, 11, 2):
    try:
        V = np.loadtxt(f'error_order_{i}.txt', delimiter=',')
        sns.distplot(np.log10(np.sqrt(V[:, column]**2)), kde_kws=dict(cumulative=True))
    except:
        pass
    
plt.show()
