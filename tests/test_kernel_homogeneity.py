"""Structure of the derivative tensor of power laws r^-a (the maths behind the
Lennard-Jones discussion in research/kernels/FINDINGS.md):

* Euler homogeneity: sum_i x_i T_a[n+e_i] = -(a+|n|) T_a[n]  (pointwise, R-dependent
  coefficients, so it reduces the *evaluation* of D but not the coefficient space);
* exponent chain: Laplacian of r^-a is a(a-1) r^-(a+2) in 3D;
* there is a constant-coefficient linear relation among the D_n (a "PDE") only
  for a = 1: the nullity of the sample matrix is 0 for every other a.
"""
import numpy as np
import pytest

from fmmgen import kernel as K
from fmmgen.kernel import multi_indices, numeric_derivative_tensor

PT = np.array([0.7, -0.4, 0.9])


@pytest.mark.parametrize("a", [1, 2, 3, 6, 12])
def test_euler_homogeneity(a):
    D = numeric_derivative_tensor(K.power_law(a), 7, PT)
    for n in multi_indices(3, 6):
        lhs = sum(PT[i] * D[tuple(e + (j == i) for j, e in enumerate(n))] for i in range(3))
        assert lhs == pytest.approx(-(a + sum(n)) * D[n], rel=1e-10, abs=1e-12)


@pytest.mark.parametrize("a", [2, 3, 6, 10])
def test_laplacian_exponent_chain(a):
    Da = numeric_derivative_tensor(K.power_law(a), 6, PT)
    Db = numeric_derivative_tensor(K.power_law(a + 2), 4, PT)
    for n in multi_indices(3, 4):
        lap = sum(Da[tuple(e + 2 * (j == i) for j, e in enumerate(n))] for i in range(3))
        assert lap == pytest.approx(a * (a - 1) * Db[n], rel=1e-9, abs=1e-12)


def _nullity(a, m, npts=400, seed=0):
    k = K.power_law(a)
    rng = np.random.default_rng(seed)
    ns = [n for n in multi_indices(3, m) if sum(n) == m]
    rows = []
    for _ in range(npts):
        R = rng.standard_normal(3)
        R /= np.linalg.norm(R)  # unit sphere: all T_n of one order have the same scaling
        D = numeric_derivative_tensor(k, m, R)
        rows.append([D[n] for n in ns])
    s = np.linalg.svd(np.array(rows), compute_uv=False)
    return int(np.sum(s < 1e-9 * s[0])), len(ns)


@pytest.mark.parametrize("m", [2, 3, 4, 5])
def test_constant_coefficient_relations_only_for_laplace(m):
    """Number of independent constant-coefficient relations among the order-m
    derivatives: m(m-1)/2 for 1/r (the harmonic reduction), none for r^-a, a != 1."""
    null1, total = _nullity(1, m)
    assert null1 == total - (2 * m + 1)
    for a in (2, 3, 6, 12):
        assert _nullity(a, m)[0] == 0, (a, m)
