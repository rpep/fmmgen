"""The kernel library (power laws, Lennard-Jones, Morse, Buckingham, Ewald
real-space, Gaussian-screened Coulomb, softened gravity): derivative tensors
against sympy and compiled far-field convergence."""
import numpy as np
import pytest

from fmmgen import kernel as K
from kernel_util import HAVE_GCC, Lib, far_field, rel_err
from test_kernel_derivatives import evaluate_D, reference

LIB = {
    "r^-1": lambda: K.power_law(1),
    "r^-3": lambda: K.power_law(3),
    "r^-6": lambda: K.power_law(6),
    "r^-12": lambda: K.power_law(12),
    "lj": lambda: K.lennard_jones(1.0, 1.0),
    "morse": lambda: K.morse(1.0, 1.5, 1.0),
    "buckingham": lambda: K.buckingham(2.0, 0.5, 1.0),
    "ewald": lambda: K.ewald_real(0.8),
    "erf": lambda: K.gaussian_screened_coulomb(0.8),
    "plummer": lambda: K.softened_power(0.2),
    "yukawa": lambda: K.yukawa(0.7),
}
PT = [0.7, -0.4, 0.9]


@pytest.mark.parametrize("name", sorted(LIB))
def test_derivative_tensor_vs_sympy(name):
    k = LIB[name]()
    assert k.method == "recurrence"
    got, ref = evaluate_D(k, 4, PT), reference(k, 4, PT)
    for n in ref:
        assert got[n] == pytest.approx(ref[n], rel=1e-8, abs=1e-9), (name, n)


@pytest.mark.skipif(not HAVE_GCC, reason="gcc not available")
@pytest.mark.parametrize("name", sorted(LIB))
def test_compiled_far_field(name, tmp_path):
    lib = Lib(LIB[name](), tmp_path, [2, 4, 6])
    errs = []
    # steep kernels (r^-6, r^-12, LJ) need well-separated, small clusters
    for p in (2, 4, 6):
        phi, F, pr, Fr = far_field(lib, p, np.random.default_rng(11), "A", radius=0.15)
        errs.append(rel_err(phi, pr))
    assert errs[2] < 0.1 * errs[0], (name, errs)
    assert errs[2] < 5e-3, (name, errs)
