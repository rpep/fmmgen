import numpy as np
import subprocess

ns = np.arange(4.0, 6.0, 0.1)
ns = np.round(10**ns)


for n in ns:
    print('n = {}'.format(n))
    p = subprocess.Popen("./lazy {} {} {} {}".format(int(n), 128, 0.5, 8), shell=True)
    stdout = p.communicate()[0]
    print('STDOUT:{}'.format(stdout))

