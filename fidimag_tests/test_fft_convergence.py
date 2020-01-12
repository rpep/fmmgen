import fidimag
import fidimag.common.constant as C
import sys

L = 12
a = 0.2715
nx, ny, nz = L, L, L
dx, dy, dz = a, a, a

J = 5.88 * C.meV
Ku = 0.41 * C.meV
mus = 3 * C.mu_B
gamma = 1.76e11
order = [4, 5, 6, 7]




theta = 0.5
ncrit = 128

mesh = fidimag.common.CuboidMesh(nx=nx, ny=ny, nz=nz, dx=dx, dy=dy, dz=dz,
                                periodicity=(False, False, False),
                                unit_length=1e-9)


sim = fidimag.atomistic.Sim(mesh, name='L_{}_FMM_convergence'.format(L), driver='llg')
sim.set_mu_s(mus)
# Add the magnetic interactions
sim.add(fidimag.atomistic.Exchange(J))
demagfft = fidimag.atomistic.Demag()
sim.add(demagfft)
sim.set_m((0, 0, 1))
sim.relax(max_steps=10000, stopping_dmdt=0.1)
sim.save_vtk()
