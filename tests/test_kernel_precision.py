"""High-order accuracy of the radial recurrence: float64 against 50-digit
arithmetic for every library kernel, and the generated array definitions for
erf(alpha r)/r, whose naive univariate derivatives cancel catastrophically."""
import pytest

from fmmgen import kernel as K
from fmmgen.kernel import numeric_derivative_tensor
from test_kernel_derivatives import evaluate_D
from test_kernel_library import LIB, PT


def _max_rel(a, b):
    scale = max(abs(v) for v in b.values())
    return max(abs(a[n] - b[n]) / scale for n in b)


@pytest.mark.parametrize("p", [12, 20])
@pytest.mark.parametrize("name", sorted(LIB))
def test_float64_recurrence_matches_50_digits(name, p):
    k = LIB[name]()
    got = numeric_derivative_tensor(k, p, PT)
    ref = numeric_derivative_tensor(k, p, PT, mp_dps=50)
    # error relative to the largest entry of each order (entries of one order are comparable)
    for m in range(1, p + 1):
        assert _max_rel(_shell(got, m), _shell(ref, m)) < 1e-10, (name, p, m)


@pytest.mark.parametrize("alpha,pt", [(0.8, PT), (0.8, [3.1, -2.2, 4.0]), (2.5, [0.1, 0.05, -0.2])])
def test_erf_derivatives_do_not_cancel(alpha, pt):
    """Regression: naive erf(a r)/r derivatives lost ~11 digits by p=20."""
    k = K.gaussian_screened_coulomb(alpha)
    got = numeric_derivative_tensor(k, 20, pt)
    ref = numeric_derivative_tensor(k, 20, pt, mp_dps=50)
    for m in range(1, 21):
        assert _max_rel(_shell(got, m), _shell(ref, m)) < 1e-10, (alpha, m)


def _shell(d, m):
    return {n: v for n, v in d.items() if sum(n) == m}


def test_erf_generated_arrays_high_order():
    k = K.gaussian_screened_coulomb(0.8)
    p = 10
    got = evaluate_D(k, p, PT)
    ref = numeric_derivative_tensor(k, p, PT, mp_dps=50)
    for m in range(1, p + 1):
        assert _max_rel(_shell(got, m), _shell(ref, m)) < 1e-12, m
