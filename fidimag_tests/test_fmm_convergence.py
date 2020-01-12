import fidimag
import fidimag.common.constant as C
import sys

L = 15
print('L = {}, L^3 = {}'.format(L, L**3))
a = 0.2715
nx, ny, nz = L, L, L
dx, dy, dz = a, a, a

J = 5.88 * C.meV
Ku = 0.41 * C.meV
mus = 3 * C.mu_B
gamma = 1.76e11
order = [4, 5, 6, 7]

theta = 0.5
print(f'theta = {theta}')

ncrit = 128

mesh = fidimag.common.CuboidMesh(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz,
                                periodicity=(False, False, False),
                                unit_length=1e-9)


sim = fidimag.atomistic.Sim(mesh, name='L_{}_FMM_convergence'.format(L), driver='llg')
sim.set_mu_s(mus)
# Add the magnetic interactions
sim.add(fidimag.atomistic.Exchange(J))
demagfmm = fidimag.atomistic.DemagFMM(order, ncrit, theta)
sim.add(demagfmm)
sim.set_m((0, 0, 1))
sim.relax()
sim.save_vtk()
