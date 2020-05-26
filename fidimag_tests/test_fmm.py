import fidimag
import fidimag.common.constant as C
import sys

L = int(sys.argv[1])
print('L = {}, L^3 = {}'.format(L, L**3))
a = 0.2715
nx, ny, nz = L, L, L
dx, dy, dz = a, a, a

J = 5.88 * C.meV
Ku = 0.41 * C.meV
mus = 3 * C.mu_B
gamma = 1.76e11
order = int(sys.argv[2])
print(f'order = {order}')
theta = float(sys.argv[3])
print(f'theta = {theta}')

ncrit = int(sys.argv[4])

mesh = fidimag.common.CuboidMesh(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz,
                                periodicity=(False, False, False),
                                unit_length=1e-9)


sim = fidimag.atomistic.Sim(mesh, name='L_{}_FMM_conf'.format(L), driver='llg')
sim.set_mu_s(mus)
# Add the magnetic interactions
sim.add(fidimag.atomistic.Exchange(J))
demagfmm = fidimag.atomistic.DemagFMM(order, ncrit, theta)
demagfft = fidimag.atomistic.Demag()
sim.add(demagfmm)
sim.add(demagfft)
sim.set_m((0, 0, 1))

import time

dt_fmm = 0
for i in range(3):
    start_fmm = time.time()
    demagfmm.compute_field()
    end_fmm = time.time()
    dt_fmm += (end_fmm - start_fmm)/3

dt_fft = 0
for i in range(3):
    start_fft = time.time()
    demagfft.compute_field()
    end_fft = time.time()
    dt_fft += (end_fft - start_fft)/3

print(f'fmm = {end_fmm - start_fmm}')
print(f'fft = {end_fft - start_fft}')
f = open(f'results/L_{L}_n_{nx*ny*nz}_order_{order}_theta_{theta}_ncrit_{ncrit}.txt', 'w')
f.write("fmm,fft\n")
f.write(f"{dt_fmm,dt_fft}")
f.close()
