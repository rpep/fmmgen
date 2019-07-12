import sys

f = open(sys.argv[1])

files = f.read().split('\n')

a = []
b = []


string = 'P2P('

for l in files:
    if 'P2P(' in l:
        lp = l.lstrip(' ')
        print(lp)
        av, bv = lp[len(string):-1].split(',')
        a.append(int(av))
        b.append(int(bv))


a, b = [list(t) for t in zip(*sorted(zip(a, b)))]

seta = set(a)
N = len(seta)




for i, j in zip(a, b):
    print('P2P({},{})'.format(i, j))


print(f'{N*(N-1)} calculations expected')
print(f'{len(a)} calculations done... ')
